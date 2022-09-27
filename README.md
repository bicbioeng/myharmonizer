# myharmonizer
Python package to harmonize and evaluate similarity between new data and an existing knowledge base in a myHarmonizer object.

## Installation
### Using Docker
```
docker run -d -e GRANT_SUDO=yes --user root --rm -p 8888:8888 -e NB_USER='jovyan' -e CHOWN_HOME=yes -w "/app" -e CHOWN_HOME_OPTS='-R' -v /home/tuyendo/USD/app/:/app --name Myharmonizer us-central1-docker.pkg.dev/nosi-usd-biofilm/nosi-usd-biofilm-arti/myharmonizer
```
Parameters:
 - -d: detach mode
 - GRANT_SUDO: grant super user for nb user
 - --user: specify docker execute user
 - --rm: automatic clean container when container stop
 - -p: export container container port for host using [host:container]
 - NB_USER: specify nb user
 - CHOWN_HOME=yes -w "/app" -e CHOWN_HOME_OPTS='-R' : Change working directory
 - -v : bind container volume to host folder
 - --name : specify container name

Requirements:
 - python <=3.9.12  

Clone git repo
```
git clone https://github.com/bicbioeng/myharmonizer
```
Install package using pip
```
pip install ./myharmonizer
```
