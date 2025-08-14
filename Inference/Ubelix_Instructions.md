# Instructions for the cluster

In general what i found easiest is to work in ipynb notebooks within vscode on the interactive ondemand option ubelix offers.

https://ondemand.hpc.unibe.ch/pun/sys/dashboard

Here you can click on interactive apps and select "(VS)Code Server". Under partition select gpu, GPU type RTX 4090 (we need all the memory we can get and the 4090 allows for the highest requestable memory). I ran with a custom instance size of 5 cores with 16GB/core. You should not have to define a CUDA version.

In the notebook run the code below to make sure you are running on a gpu capable core. Sometimes gpu nodes have issues so be sure to check you are actually have acess to it.
'import torch
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print("Using:", device)'


