# Welcome to the Jan 2024 ESS McStas Workshop

## Workshop programme:
[![Workshop programme](pics/programme.png)](https://docs.google.com/spreadsheets/d/1x_YmHhAyZquYxcbj0iaY-HADBu12bvtAuzNeIsgVlHA/edit?usp=sharing)

## Prerequisites, 
- Local installation of [docker](https://www.docker.com/products/docker-desktop)
- Account on ESS DMSC server
- (Optionally [XQuartz](https://www.xquartz.org) if on macOS or [Xming](https://sourceforge.net/projects/xming/files/latest/download) if on Windows)

## Starting the docker
- ```docker run -p 8888:8888 docker.io/mccode/mcstas-2.7.1-3.1-scipp:1.0``` 
- (add ```-v /some/folder:/home/jovyan/otherfolder``` before the image name to map a folder from your local machine to the docker image)
- Connect to the URL communicated by the docker (will open a JupyterLab that includes a link to a browser-based X-session)
![Jupyter started](pics/Jupyterlab.png)
- Pressing the Desktop link will open an XFCE session
![Finding McStas](pics/XFCE_McStas_launchers.png)
![McStas started](pics/McStas_3.1_started.png)

## Zoom access
A zoom link will be shared later
