# MicrobiomeHub

## miniconda
https://github.com/luozhy88/MicrobiomeHub/blob/main/install/01_install_miniconda.sh

## microbiome.net(mamba)
https://github.com/luozhy88/MicrobiomeHub/blob/main/install/02_install_mamba.sh

## R package
https://github.com/luozhy88/MicrobiomeHub/blob/main/install/03_install_r_packages.r

## shiny-server 

### 重新定义R版本
/etc/systemd/system/shiny-server.service中添加：  
Environment="R=/home/test1/miniconda3/envs/microbiome.net/bin/R"  #30site  
sudo systemctl stop shiny-server  
sudo systemctl daemon-reload   
sudo systemctl start shiny-server  

### 更改工作目录
sudo vim /etc/shiny-server/shiny-server.conf   
sudo systemctl restart shiny-server   
sudo systemctl restart shiny-server.service   


sudo systemctl stop shiny-server  
sudo systemctl daemon-reload   
sudo systemctl start shiny-server  
sudo systemctl restart shiny-server   
sudo systemctl restart shiny-server.service   



## Note
log:/var/log/shiny-server/*.log  
test shiny:https://github.com/luozhy88/MicrobiomeHub/blob/main/hello/server.R
