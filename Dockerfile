# Base image https://hub.docker.com/u/rocker/
FROM rocker/shiny:4.3.1

## update system libraries and install wget
RUN sudo apt-get update -y
RUN sudo apt-get install -y wget
RUN sudo apt-get install -y make
RUN sudo apt-get install -y build-essential
RUN sudo apt-get install -y gcc-multilib
RUN sudo apt-get install -y systemctl

# downlad, unpack, install openBugs
RUN wget -P /tmp https://web.archive.org/web/20230611171942/https://www.mrc-bsu.cam.ac.uk/wp-content/uploads/2018/04/OpenBUGS-3.2.3.tar.gz
RUN tar zxvf tmp/OpenBUGS-3.2.3.tar.gz
WORKDIR /OpenBUGS-3.2.3
RUN ./configure
RUN make
RUN make install
WORKDIR /

# java files
RUN sudo apt-get install -y default-jre
RUN sudo apt-get install -y default-jdk
RUN sudo R CMD javareconf
#RUN sudo apt install -y r-cran-rjava

# install latex
RUN apt-get install texlive-latex-recommended texlive-fonts-recommended -y
RUN apt-get install texlive-latex-extra -y

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
RUN Rscript -e 'install.packages("rmarkdown")'
RUN Rscript -e 'install.packages("tinytex")'
RUN Rscript -e 'tinytex::install_tinytex()'
RUN Rscript -e 'install.packages("knitr")'

# delete default apps
WORKDIR /srv/shiny-server
RUN rm -r 01_hello 02_text 03_reactivity 04_mpg 05_sliders 06_tabsets 07_widgets 08_html 09_upload 10_download 11_timer
RUN chmod -R 755 /srv/shiny-server

# expose port
EXPOSE 3838

# CMD R -e 'shiny::runApp(port = 3838, host = "0.0.0.0")'
# app should run from 'inherited' RUN command from rocker/shiny
