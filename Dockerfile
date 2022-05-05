# Base image https://hub.docker.com/u/rocker/
FROM rocker/shiny:4.1.0

## update system libraries and install wget
RUN sudo apt update
RUN sudo apt install -y wget
RUN sudo apt install -y make
RUN sudo apt install -y build-essential
RUN sudo apt install -y gcc-multilib
RUN sudo apt install -y systemctl

# downlad, unpack, install openBugs
RUN wget -P /tmp https://www.mrc-bsu.cam.ac.uk/wp-content/uploads/2018/04/OpenBUGS-3.2.3.tar.gz
RUN tar zxvf tmp/OpenBUGS-3.2.3.tar.gz
WORKDIR /OpenBUGS-3.2.3
RUN ./configure
RUN make
RUN make install
WORKDIR /

# java files
RUN sudo apt update
RUN sudo apt install -y default-jre
RUN sudo apt install -y default-jdk
RUN sudo R CMD javareconf
#RUN sudo apt install -y r-cran-rjava


# copy necessary files
COPY ./ /srv/shiny-server/

# install packages
RUN Rscript -e 'install.packages("shinycssloaders")'
RUN Rscript -e 'install.packages("shinyBS")'
RUN Rscript -e 'install.packages("Matrix")'
RUN Rscript -e 'install.packages("magic")'
RUN Rscript -e 'install.packages("R2OpenBUGS")'
RUN Rscript -e 'install.packages("rJava")'
RUN Rscript -e 'install.packages("xlsx")'

# expose port
EXPOSE 3838

# CMD R -e 'shiny::runApp(port = 3838, host = "0.0.0.0")'
# app should run from 'inherited' RUN command from rocker/shiny