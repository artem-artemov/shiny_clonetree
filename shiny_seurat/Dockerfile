FROM rocker/shiny

RUN apt-get update && apt-get install --yes \    
       libssl-dev \
       libcurl4-openssl-dev \
       libxml2-dev \
       libpng-dev 
#    R -e "install.packages(c('class', 'Seurat'))"; \
#    apt-get clean;

RUN R -e "install.packages(c('Seurat', 'ggplot2', 'scales', 'visNetwork'),dependencies=TRUE, repos='http://cran.rstudio.com/')"

