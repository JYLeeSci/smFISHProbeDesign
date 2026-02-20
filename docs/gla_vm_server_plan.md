# Glasgow VM Server — Deployment Plan

Deploy the smFISHProbeDesign Streamlit app on an Oracle Linux 9.6 university VM.

---

## VM Specifications

| Resource | Value |
|----------|-------|
| OS | Oracle Linux Server 9.6 (RHEL-family, `el9`) |
| CPU | Intel Xeon Silver 4310 @ 2.10 GHz — 2 cores / 4 threads |
| RAM | 3.4 GB total |
| Swap | 3.9 GB |
| Root disk (`/`) | 70 GB (61 GB free) |
| Home disk (`/home`) | 425 GB (421 GB free) |

## Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Package manager | micromamba | Fast, lightweight, no Anaconda licensing |
| Install location | `~/GitHub/smFISHProbeDesign` | User's home directory, easy `git pull` updates |
| Process manager | systemd | Auto-restart, survives reboots, journal logging |
| Network exposure | Direct Streamlit on port 8501 | University-internal only, simpler than nginx |
| RepeatMasker | Skip for now | Users can upload pre-masked files; add later if needed |
| Swap | Keep existing 3.9 GB | Already sufficient — genome masking peaks ~3.5 GB |

## Disk Space Budget

| Component | Size | Location |
|-----------|------|----------|
| Repository (code + pseudogene FASTAs) | ~200 MB | `~/GitHub/smFISHProbeDesign` |
| Conda environment (Python, Bowtie, Streamlit) | ~500 MB | `~/micromamba/envs/probedesign` |
| Pseudogene bowtie indices (3 species) | ~85 MB | `~/GitHub/smFISHProbeDesign/bowtie_indexes` |
| Genome indices (human + mouse + drosophila) | ~6.7 GB | `~/GitHub/smFISHProbeDesign/bowtie_indexes` |
| **Total** | **~7.5 GB** | All under `/home` (421 GB free) |

---

## Step-by-step Deployment

### Step 1 — SSH into the VM

```bash
ssh <your-username>@<vm-ip-address>
```

Verify the OS and resources:

```bash
cat /etc/os-release          # Should show Oracle Linux 9.6
free -h                      # Confirm ~3.4 GB RAM + ~3.9 GB swap
df -h /home                  # Confirm ~421 GB free
```

---

### Step 2 — Install System Prerequisites

Install compilers and tools needed by conda packages and bowtie index building.

```bash
sudo dnf install -y git tar unzip wget curl gcc gcc-c++ make bzip2
```

Verify git:

```bash
git --version
```

---

### Step 3 — Install micromamba

micromamba is a fast, standalone conda package manager. Install it into your home directory (no root required for this step).

```bash
# Download and install micromamba
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
```

During install, accept the defaults:
- Install location: `~/micromamba` (or `~/.local/bin/micromamba`)
- Shell init: Yes (adds to `~/.bashrc`)

After install, reload your shell:

```bash
source ~/.bashrc
micromamba --version
```

Expected output: `2.x.x` or similar.

**Configure channels** (one-time):

```bash
micromamba config append channels conda-forge
micromamba config append channels bioconda
micromamba config set channel_priority strict
```

---

### Step 4 — Clone the Repository

```bash
mkdir -p ~/GitHub
cd ~/GitHub
git clone https://github.com/JYLeeSci/smFISHProbeDesign.git
cd smFISHProbeDesign
```

---

### Step 5 — Run the Setup Script Inside a screen Session

The setup script downloads ~6.7 GB of genome indices and takes 15–30+ minutes depending on network speed. **Run it inside `screen`** so the process survives SSH disconnection. If your connection drops mid-download, you can reconnect and reattach to find it still running.

**5a. Install `screen`** (if not already present):

```bash
sudo dnf install -y screen
```

**5b. Start a named screen session:**

```bash
screen -S probedesign_setup
```

You are now inside a persistent session. The terminal looks the same — the difference is the session lives on the server independently of your SSH connection.

**5c. Run the setup script:**

```bash
cd ~/GitHub/smFISHProbeDesign
chmod +x setup_all.sh
./setup_all.sh
```

**Essential `screen` commands:**

| Action | Keys |
|--------|------|
| Detach (leave running) | `Ctrl-A` then `D` |
| Reattach after reconnecting | `screen -r probedesign_setup` |
| List all sessions | `screen -ls` |
| Kill session when done | `exit` (inside session) |

If your SSH connection drops, simply SSH back in and run `screen -r probedesign_setup` to reattach.

**What the script does:**

1. **Create `probedesign` conda environment** (~500 MB download)
   - Python 3.11, Bowtie 1.3.1, Click, Streamlit, Pandas
2. **Install the `probedesign` package** (editable mode)
3. **Build pseudogene indices** from shipped FASTAs (~1 min)
   - `humanPseudo`, `mousePseudo`, `drosophilaPseudo`
4. **Download genome indices** from AWS (~6.7 GB, ~12–30 min depending on network)
   - Human GRCh38 (~3.5 GB), Mouse mm10 (~3.1 GB), Drosophila BDGP6 (~176 MB)
5. **Run validation tests** (7 checks)

**If you only need pseudogene masking** (skip the 6.7 GB genome download):

```bash
./setup_all.sh --skip-genome
```

You can always re-run `./setup_all.sh` later — it skips already-completed steps, so it is safe to resume an interrupted run.

**Expected output** (all tests pass):

```
Tests: 7/7 passed
All checks passed!
```

When done, exit the screen session:

```bash
exit
```

---

### Step 6 — Verify the Installation

Activate the environment and run quick checks.

```bash
micromamba activate probedesign

# Check tools are available
python --version           # Python 3.11.x
bowtie --version           # bowtie-align-s version 1.3.1
streamlit --version        # 1.x.x
probedesign --version      # 0.1.0
```

**Verify bowtie index auto-detection:**

Both the CLI and the Streamlit app derive the index directory automatically from the script's own location using `Path(__file__).resolve()`, not from the current working directory. This means the path is always `~/GitHub/smFISHProbeDesign/bowtie_indexes` regardless of where you `cd` before running.

Confirm the indices are where the app expects them:

```bash
ls ~/GitHub/smFISHProbeDesign/bowtie_indexes/*.1.ebwt 2>/dev/null | head -5
# Should list: humanPseudo.1.ebwt, mousePseudo.1.ebwt, drosophilaPseudo.1.ebwt

ls ~/GitHub/smFISHProbeDesign/bowtie_indexes/*.1.bt2 2>/dev/null | head -5
# Should list: GCA_000001405.15_GRCh38_no_alt_analysis_set.1.bt2, mm10.1.bt2, drosophila.1.bt2
```

Test CLI probe design with pseudogene masking to confirm the index is found:

```bash
cd ~/GitHub/smFISHProbeDesign

# No masking — fast sanity check
probedesign design test_cases/CDKN1A_32/CDKN1A.fa --probes 32 \
  --repeatmask-file test_cases/CDKN1A_32/CDKN1A_repeatmasked.fa \
  -o /tmp/test_cdkn1a --quiet
cat /tmp/test_cdkn1a_oligos.txt | wc -l   # Should print: 32

# With pseudogene masking (verifies index auto-detection)
probedesign design test_cases/KRT19_withUTRs/KRT19_withUTRs.fa --probes 32 \
  --pseudogene-mask --quiet -o /tmp/test_krt19
echo "Exit code: $?"   # Should print: Exit code: 0
```

In the Streamlit UI, the "Bowtie index directory" sidebar field is pre-filled with the same auto-detected path. The `check_prerequisites()` function validates the index files exist before any design run and shows a warning banner if anything is missing — it will never silently fail.

**Quick Streamlit smoke test** (verify it starts, then Ctrl+C):

```bash
streamlit run streamlit_app/app.py --server.port 8501 --server.address 0.0.0.0
# You should see: "You can now view your Streamlit app in your browser."
# Press Ctrl+C to stop — the systemd service will manage it from here
```

---

### Step 7 — Create a systemd Service

A systemd unit ensures the app starts on boot, restarts on crash, and logs to journald.

**7a. Find the full path to Streamlit in the conda environment:**

```bash
# While probedesign env is activated:
which streamlit
# Example output: /home/<username>/micromamba/envs/probedesign/bin/streamlit
```

Note this path — you will use it in the service file below.

**7b. Create the service file:**

```bash
sudo tee /etc/systemd/system/probedesign.service > /dev/null << 'SERVICEEOF'
[Unit]
Description=ProbeDesign Streamlit Web App
After=network.target

[Service]
Type=simple
User=<your-username>
Group=<your-username>

# Working directory — the cloned repo
WorkingDirectory=/home/<your-username>/GitHub/smFISHProbeDesign

# Full path to streamlit inside the conda environment
# IMPORTANT: Use the output from `which streamlit` above
ExecStart=/home/<your-username>/micromamba/envs/probedesign/bin/streamlit run \
    streamlit_app/app.py \
    --server.port 8501 \
    --server.address 0.0.0.0 \
    --server.headless true \
    --browser.gatherUsageStats false

# Restart policy
Restart=on-failure
RestartSec=5

# Environment — ensure the conda env's bin is in PATH so bowtie is found
Environment="PATH=/home/<your-username>/micromamba/envs/probedesign/bin:/usr/local/bin:/usr/bin:/bin"
Environment="HOME=/home/<your-username>"

# Resource limits
LimitNOFILE=65536

# Logging
StandardOutput=journal
StandardError=journal
SyslogIdentifier=probedesign

[Install]
WantedBy=multi-user.target
SERVICEEOF
```

**IMPORTANT:** Replace every `<your-username>` with your actual username. There are 6 occurrences.

**7c. Enable and start the service:**

```bash
sudo systemctl daemon-reload
sudo systemctl enable probedesign        # Start on boot
sudo systemctl start probedesign         # Start now
```

**7d. Check the service is running:**

```bash
sudo systemctl status probedesign
```

Expected output should show `active (running)`. If it fails, check logs:

```bash
sudo journalctl -u probedesign -f        # Follow live logs
sudo journalctl -u probedesign --no-pager -n 50   # Last 50 lines
```

---

### Step 8 — Configure the Firewall

Oracle Linux 9 uses `firewalld`. Open port 8501 for the Streamlit app:

```bash
sudo firewall-cmd --permanent --add-port=8501/tcp
sudo firewall-cmd --reload

# Verify
sudo firewall-cmd --list-ports
# Should include: 8501/tcp
```

If `firewalld` is not running (check with `sudo systemctl status firewalld`), you may need to start it or use `iptables` instead:

```bash
# Only if firewalld is not available:
sudo iptables -A INPUT -p tcp --dport 8501 -j ACCEPT
sudo iptables-save | sudo tee /etc/sysconfig/iptables
```

**University network security group**: Check with your IT department whether the VM's network security group (or equivalent) allows inbound TCP connections on port 8501. This is separate from the OS firewall. If port 8501 is blocked at the network level, the firewall commands above are not sufficient.

---

### Step 9 — Verify the Deployment

From another machine on the university network, open a browser and navigate to:

```
http://<vm-ip-address>:8501
```

You should see the ProbeDesign Streamlit interface.

**Test a probe design run:**

1. Paste or upload a FASTA sequence (e.g., from `test_cases/CDKN1A_32/CDKN1A.fa`)
2. Set probes to 32
3. Enable Pseudogene mask and Genome mask
4. Click **Design Probes**
5. Verify results appear and files are downloadable

**Monitor memory during genome masking** (in a separate SSH session):

```bash
watch -n 1 free -m
```

During genome masking, expect RAM + swap to be heavily used (~3.5 GB for human genome bowtie queries). The existing 3.9 GB swap should handle this, but observe to confirm the process does not get OOM-killed.

---

## Maintenance & Updates

### Updating the App (git pull)

When new features are pushed to the repository:

```bash
# SSH into the VM
cd ~/GitHub/smFISHProbeDesign

# Pull latest changes
git pull origin master

# If environment.yml changed (new dependencies), update the env:
micromamba activate probedesign
micromamba env update -n probedesign -f environment.yml --prune -y

# If pyproject.toml changed, reinstall the package:
pip install -e .

# Restart the service to pick up changes
sudo systemctl restart probedesign

# Verify
sudo systemctl status probedesign
```

**Quick update** (if only Python code changed, no new dependencies):

```bash
cd ~/GitHub/smFISHProbeDesign
git pull origin master
sudo systemctl restart probedesign
```

### Viewing Logs

```bash
# Live log stream
sudo journalctl -u probedesign -f

# Last 100 lines
sudo journalctl -u probedesign --no-pager -n 100

# Logs from today only
sudo journalctl -u probedesign --since today

# Search for errors
sudo journalctl -u probedesign --no-pager | grep -i error
```

### Restarting the Service

```bash
sudo systemctl restart probedesign      # Restart
sudo systemctl stop probedesign         # Stop
sudo systemctl start probedesign        # Start
sudo systemctl status probedesign       # Check status
```

### Monitoring Memory

```bash
# One-shot check
free -m

# Continuous monitoring (update every 2 seconds)
watch -n 2 free -m

# Check if OOM killer has acted
sudo dmesg | grep -i "out of memory"
sudo journalctl --since today | grep -i "oom"
```

---

## Troubleshooting

### Service fails to start

```bash
# Check the logs
sudo journalctl -u probedesign --no-pager -n 50

# Common issues:
# 1. Wrong path to streamlit → fix ExecStart in service file
# 2. Wrong username → fix User/Group and paths in service file
# 3. Missing bowtie → check PATH environment in service file
```

After editing the service file:

```bash
sudo systemctl daemon-reload
sudo systemctl restart probedesign
```

### "bowtie: command not found" in service

The systemd service runs without shell initialization, so conda's PATH isn't set automatically. The `Environment="PATH=..."` line in the service file must include the conda env's bin directory. Verify:

```bash
# Find where bowtie is installed
micromamba activate probedesign
which bowtie
# Example: /home/<username>/micromamba/envs/probedesign/bin/bowtie

# The PATH in the service file must include that directory:
# Environment="PATH=/home/<username>/micromamba/envs/probedesign/bin:..."
```

### App accessible on VM but not from other machines

1. **OS firewall**: Verify port 8501 is open:
   ```bash
   sudo firewall-cmd --list-ports
   ```
2. **Network security group**: Contact IT — the VM's cloud/network security group may block the port.
3. **Streamlit binding**: Ensure `--server.address 0.0.0.0` (not `127.0.0.1`) in the service file.

### Genome masking causes OOM kill

Symptoms: The app hangs or crashes during genome masking. Check:

```bash
sudo dmesg | tail -20   # Look for "Out of memory: Killed process"
```

If this happens, the 3.4 GB RAM + 3.9 GB swap may not be enough under heavy load. Options:

1. **Add more swap** (easiest):
   ```bash
   sudo fallocate -l 4G /swapfile
   sudo chmod 600 /swapfile
   sudo mkswap /swapfile
   sudo swapon /swapfile

   # Make persistent across reboots
   echo '/swapfile none swap sw 0 0' | sudo tee -a /etc/fstab

   # Verify
   free -h   # Should now show ~7.9 GB swap total
   ```

2. **Limit concurrent users**: Only one genome masking job should run at a time. Streamlit's default behavior (one Python process per browser tab) makes this unlikely unless multiple users run genome masking simultaneously.

3. **Request more RAM from IT** if genome masking is a critical feature.

### Streamlit shows "Please wait..." forever

```bash
# Check if the process is running
sudo systemctl status probedesign

# Check if port is listening
ss -tlnp | grep 8501

# Try accessing locally on the VM
curl -s http://localhost:8501 | head -5
```

### Updating fails (git pull conflicts)

If someone edited files directly on the VM:

```bash
cd ~/GitHub/smFISHProbeDesign

# Option 1: Stash local changes, pull, then reapply
git stash
git pull origin master
git stash pop

# Option 2: Discard local changes entirely (destructive)
git checkout -- .
git pull origin master
```

---

## Optional: Install RepeatMasker (Later)

If you later decide to enable automatic repeat masking, follow these steps. This requires ~65 GB free disk space during installation (~56 GB permanent).

### 1. Install RepeatMasker via conda

```bash
micromamba activate probedesign
micromamba install -c bioconda -c conda-forge repeatmasker -y
```

### 2. Download the Dfam database (Mammalia partition)

This is the large download — ~8.9 GB compressed, ~56 GB extracted.

```bash
# Find the RepeatMasker Libraries directory
RMPATH=$(dirname $(which RepeatMasker))/../share/RepeatMasker/Libraries/famdb
echo "RepeatMasker famdb dir: $RMPATH"

cd "$RMPATH"

# Download Partition 0 (root, always needed) — if not already present
# It's usually included with the conda install, but verify:
ls dfam39_full.0.h5 2>/dev/null || {
    curl -O https://www.dfam.org/releases/current/families/FamDB/dfam39_full.0.h5.gz
    gunzip dfam39_full.0.h5.gz
}

# Download Partition 7 (Mammalia — needed for human, mouse, rat)
curl -O https://www.dfam.org/releases/current/families/FamDB/dfam39_full.7.h5.gz
gunzip dfam39_full.7.h5.gz    # This creates a ~56 GB file
```

### 3. Verify RepeatMasker

```bash
micromamba activate probedesign
RepeatMasker -v       # Should print version
RepeatMasker -species human -dir /tmp/rm_test test_cases/CDKN1A_32/CDKN1A.fa
ls /tmp/rm_test/      # Should contain CDKN1A.fa.masked
```

### 4. Restart the service

```bash
sudo systemctl restart probedesign
```

RepeatMasker will now be detected automatically by the app's prerequisite check. The "Auto RepeatMasker" option will become available in the sidebar.

---

## Optional: Nginx Reverse Proxy (Later)

If you later want to add HTTPS, a clean URL (no `:8501` port), or basic authentication, you can put Nginx in front of Streamlit.

### 1. Install Nginx

```bash
sudo dnf install -y nginx
sudo systemctl enable --now nginx
```

### 2. Create Nginx config

```bash
sudo tee /etc/nginx/conf.d/probedesign.conf > /dev/null << 'NGINXEOF'
server {
    listen 80;
    server_name <vm-hostname-or-ip>;

    location / {
        proxy_pass http://127.0.0.1:8501;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
        proxy_read_timeout 86400;
    }
}
NGINXEOF
```

### 3. Update Streamlit to bind to localhost only

Edit the systemd service file to change `--server.address 0.0.0.0` to `--server.address 127.0.0.1` so only Nginx can reach Streamlit.

### 4. Open port 80 and close 8501

```bash
sudo firewall-cmd --permanent --add-service=http
sudo firewall-cmd --permanent --remove-port=8501/tcp
sudo firewall-cmd --reload
```

### 5. Restart services

```bash
sudo systemctl restart probedesign
sudo systemctl restart nginx
```

Access via `http://<vm-hostname-or-ip>` (no port number).

---

## Quick Reference

| Task | Command |
|------|---------|
| SSH in | `ssh <username>@<vm-ip>` |
| Start screen session | `screen -S probedesign_setup` |
| Detach from screen | `Ctrl-A` then `D` |
| Reattach to screen | `screen -r probedesign_setup` |
| Check service | `sudo systemctl status probedesign` |
| Restart service | `sudo systemctl restart probedesign` |
| View live logs | `sudo journalctl -u probedesign -f` |
| Update app | `cd ~/GitHub/smFISHProbeDesign && git pull origin master && sudo systemctl restart probedesign` |
| Update env | `micromamba activate probedesign && micromamba env update -n probedesign -f environment.yml --prune -y && pip install -e . && sudo systemctl restart probedesign` |
| Check memory | `free -m` |
| Check port | `ss -tlnp \| grep 8501` |
| Check firewall | `sudo firewall-cmd --list-ports` |
| Verify index dir | `ls ~/GitHub/smFISHProbeDesign/bowtie_indexes/*.1.bt2 \| wc -l` |
| App URL | `http://<vm-ip>:8501` |
