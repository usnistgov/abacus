#!/bin/sh

set -e

rm master.zip
rm -r abacus-main

wget https://github.com/usnistgov/abacus/archive/refs/heads/main.zip
unzip master.zip

sudo docker build -t abacus ./abacus-main
sudo docker rm -f $(sudo docker ps -a -q)
sudo docker image prune -f
sudo docker run -d -p 80:3838 abacus
