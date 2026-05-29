"""
Hybrid RSA+AES-GCM encryption utilities for validator-miner communication.

The validator generates a fresh RSA keypair per miner query:
  - public key  (encryption_key) → sent to miner inside the synapse
  - private key (decryption_key) → kept by the validator, never shared

The miner encrypts its response with the public key using hybrid encryption:
  - a random 256-bit AES-GCM key encrypts the plaintext
  - the AES key is then encrypted with the RSA public key (OAEP/SHA-256)

Because each miner receives a unique public key, no miner can re-use or
re-submit another miner's encrypted payload.
"""

import base64
import json
import os

from cryptography.hazmat.primitives import hashes, serialization
from cryptography.hazmat.primitives.asymmetric import padding, rsa
from cryptography.hazmat.primitives.ciphers.aead import AESGCM


def generate_keypair() -> tuple[str, str]:
    """
    Generate a 2048-bit RSA keypair.

    Returns:
        (public_key_pem, private_key_pem) as PEM strings.
    """
    private_key = rsa.generate_private_key(public_exponent=65537, key_size=2048)
    public_pem = private_key.public_key().public_bytes(
        serialization.Encoding.PEM,
        serialization.PublicFormat.SubjectPublicKeyInfo,
    ).decode()
    private_pem = private_key.private_bytes(
        serialization.Encoding.PEM,
        serialization.PrivateFormat.TraditionalOpenSSL,
        serialization.NoEncryption(),
    ).decode()
    return public_pem, private_pem


def encrypt(public_key_pem: str, plaintext: str) -> str:
    """
    Encrypt *plaintext* with the RSA public key using hybrid RSA+AES-GCM.

    Returns a JSON string with base64-encoded fields:
      {
        "key":   <base64 RSA-OAEP-encrypted AES key>,
        "nonce": <base64 AES-GCM nonce (12 bytes)>,
        "data":  <base64 AES-GCM ciphertext + authentication tag>
      }
    """
    public_key = serialization.load_pem_public_key(public_key_pem.encode())

    aes_key = os.urandom(32)   # 256-bit AES key
    nonce = os.urandom(12)     # 96-bit GCM nonce

    ciphertext = AESGCM(aes_key).encrypt(nonce, plaintext.encode(), None)

    encrypted_key = public_key.encrypt(
        aes_key,
        padding.OAEP(
            mgf=padding.MGF1(algorithm=hashes.SHA256()),
            algorithm=hashes.SHA256(),
            label=None,
        ),
    )

    return json.dumps({
        "key":   base64.b64encode(encrypted_key).decode(),
        "nonce": base64.b64encode(nonce).decode(),
        "data":  base64.b64encode(ciphertext).decode(),
    })


def decrypt(private_key_pem: str, encrypted_payload: str) -> str:
    """
    Decrypt a payload produced by :func:`encrypt`.

    Args:
        private_key_pem:    PEM-encoded RSA private key.
        encrypted_payload:  JSON string returned by :func:`encrypt`.

    Returns:
        The original plaintext string.
    """
    private_key = serialization.load_pem_private_key(
        private_key_pem.encode(), password=None
    )

    payload = json.loads(encrypted_payload)
    encrypted_key = base64.b64decode(payload["key"])
    nonce = base64.b64decode(payload["nonce"])
    ciphertext = base64.b64decode(payload["data"])

    aes_key = private_key.decrypt(
        encrypted_key,
        padding.OAEP(
            mgf=padding.MGF1(algorithm=hashes.SHA256()),
            algorithm=hashes.SHA256(),
            label=None,
        ),
    )

    return AESGCM(aes_key).decrypt(nonce, ciphertext, None).decode()
