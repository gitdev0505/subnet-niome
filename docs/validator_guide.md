## Validator Setup Guide

To contribute as a Validator in the NIOME(SN55), you must prepare your machine, install the required packages, and register your identity on Bittensor. Validators play a critical role in evaluating miner outputs and maintaining subnet integrity.

### 1. Prerequisites

Before starting, ensure your system meets the requirements and has the core dependencies installed.

* **Operating System:** Ubuntu 22.04 or similar Linux distribution is generally recommended for optimal compatibility. Mining is not supported on Windows.
* **Python:** Python 3.12
* **Git:** required for cloning the repository
* **Hardware:** 
   - vCPU :  32
   - GPU : unnecessary
   - Memory : 64GB recommended
   - Storage : 1TB SSD recommended
   - 3rd Party API : unnecessary
   - Port Forwarding : standard

Validaors must maintain high uptime and stable networking, as they responsible for scoring miners and producing consensus signals

### 2. Environment Setup

This section walks you through cloning the subnet-niome repository and installing the required packages.

1. **Clone the Repository:**
   **Bash**

   ```
   git clone https://github.com/genomesio/subnet-niome.git
   cd subnet-niome
   ```
2. **Create a Virtual Environment (Recommended):**
   **Bash**

   ```
   python3 -m venv venv
   source venv/bin/activate
   ```
3. **Install Dependencies:** Install the required Python packages and register the local package for execution.
   **Bash**

   ```
   python3 -m pip install -r requirements.txt
   ```
4. **Install Docker:** Install the Docker for PharmCAT
   **Bash**

   ***Set up the repository***
   ```
   sudo apt-get update 
   sudo apt-get install -y ca-certificates curl gnupg
   sudo install -m 0755 -d /etc/apt/keyrings
   curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /etc/apt/keyrings/docker.gpg
   sudo chmod a+r /etc/apt/keyrings/docker.gpg 
   ```

   ***Add Docker's APT source***
   ```
   echo \ 
   "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.gpg] https://download.docker.com/linux/ubuntu \ 
   $(. /etc/os-release && echo "$VERSION_CODENAME") stable" | \
   sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
   ```
   
   ***Install Docker Packages***
   ```
   sudo apt-get update 
   sudo apt-get install -y docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin
   ```

   ***Verify it works***
   ```
   docker run --rm hello-world
   ```

### 3. Running the Validator

Once your hotkey is registered, you can start your Validator. 

1. **Run the Validator Script:** The core command to launch a validator neuron requires specifying your wallet and hotkey names, the network, and the subnet ID (`--netuid 55`).
   
   In your current subnet-niome project path

   **Bash**
   ```
   export PYTHONPATH="$PYTHONPATH:$(pwd)                                                               
   ```

   ```
   python neurons/validator.py \
   --netuid 55 \
   --subtensor.network finney \
   --wallet.name your_coldkey \
   --wallet.hotkey your_hotkey \
   --wandb.api_key your_api_key \ 
   --logging.debug
   ```

3. **Keep it Running:** Use a process manager like **`pm2`** or **`tmux`** to ensure your validator stays online. Validators must remain responsive to score miners and participate in consensus; downtime reduces emissions