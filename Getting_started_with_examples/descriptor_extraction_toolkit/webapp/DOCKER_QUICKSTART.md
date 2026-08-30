# DescriPytor Docker Quickstart

This folder contains a Docker wrapper for the descriptor-extraction webapp.

Docker gives the app a reproducible Linux chemistry environment with Streamlit,
RDKit, ASE, Morfeus, AQME, Mordred, MolTop, xTB, and the local DescriPyTor
source code.

## 1. Install WSL2 And Docker

On Windows, Docker Desktop uses WSL2 for Linux containers.

1. Open **PowerShell as Administrator**.
2. Install WSL2:

   ```powershell
   wsl --install
   ```

3. Reboot if Windows asks.
4. If Ubuntu opens, create the Linux user, then type:

   ```bash
   exit
   ```

5. Back in Windows PowerShell, confirm WSL is installed:

   ```powershell
   wsl -l -v
   ```

   Ubuntu should appear with `VERSION 2`.

6. Download and install Docker Desktop from
   https://www.docker.com/products/docker-desktop/
7. Start Docker Desktop.
8. Enable Ubuntu integration:

   ```text
   Docker Desktop -> Settings -> Resources -> WSL Integration
   ```

9. Wait until Docker Desktop says the engine is running.
10. Open a new PowerShell terminal.

Check it:

```powershell
docker version
```

The command should show both `Client` and `Server`.

## 2. Build and run from the repo

`docker-compose.yml` lives at the **repository root**, not in this `webapp`
folder. From the clone:

```powershell
cd path\to\DescriPyTor_to_upload
docker compose up --build
```

Open:

```text
http://localhost:8503
http://localhost:7432/visual
```

The compose file maps host port `8503` to container port `8501` by default to
avoid conflicts with other Streamlit apps. To use a different host port:

```powershell
$env:DESCRIPYTOR_HOST_PORT=8504
docker compose up --build
```

Then open `http://localhost:8504`.

## 3. Use the app

1. Go to **Data** and upload `.feather` and/or `.xyz` files.
2. Optionally use **Pick atoms** to inspect a molecule and export atom choices.
3. Go to **Configure & run**.
4. Choose engines or upload/edit a config.
5. Download the generated `run_config.json` if you want to keep the exact run.
6. Click **Run extraction**.
7. Go to **Results** and download `merged_features.csv`.

## 4. Send it to someone else

There are two practical ways.

### Option A: send the repo

Send the repository or push it to GitHub. The other person installs
Docker Desktop and runs:

```powershell
cd path\to\DescriPyTor_to_upload
docker compose up --build
```

This is best if they may edit the code.

### Option B: send a built image

Build once:

```powershell
docker compose build
```

Save the image:

```powershell
docker save descripytor-webapp:latest -o descripytor-webapp.tar
```

Send `descripytor-webapp.tar`.

On the other computer:

```powershell
docker load -i descripytor-webapp.tar
docker run --rm -p 8503:8501 descripytor-webapp:latest
```

Open:

```text
http://localhost:8503
```

This image contains the app source. The user can upload input files through the
browser and download the final CSV.

## 5. Share inside the lab

Run it on one lab machine:

```powershell
docker compose up --build
```

Find that machine's local IP address, then labmates on the same network open:

```text
http://MACHINE-IP:8503
```

Keep this on a trusted network. The app currently has no login screen.
