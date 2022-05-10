# NIST Abacus Application

This R Shiny application is the result of a collaboration between CSD and SED. Abacus automates experimental design and statistical analysis for common procedures used in chemical analyis. 

## Using Abacus Locally

The easiest and most consistent way to run and use Abacus locally is by first installing Docker (https://docs.docker.com/get-docker/).

Once Docker is installed, download the Abacus repository (at github.com/usnistgov/abacus, click 'Code' and then 'Download ZIP'), and navigate your terminal to the main directory of the project (the same level as Dockerfile). Then, run the following command to build the image:
```
docker build -t abacus .
```
To run the container, run the following command:
```
docker run -d -p 8080:3838 --name my_container abacus
```
(-d for 'detached', -p specifies the port mapping, '--name' gives a name to the running container, and 'abacus' tells docker which image to build the container from.) Then the app should be visible at localhost:8080 (accessed via your web browser).

To stop and remove the running container, run the following:
```
docker rm -f my_container
```

Alternatively, you can run, stop, and remove the container using the Docker Desktop user interface.

For questions or bug reports regarding the software, please contact david.newton@nist.gov.
