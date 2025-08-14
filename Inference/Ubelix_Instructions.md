# Instructions for the cluster

In general what i found easiest is to work in ipynb notebooks within vscode on the interactive ondemand option ubelix offers.

https://ondemand.hpc.unibe.ch/pun/sys/dashboard

Here you can click on interactive apps and select "(VS)Code Server". Under partition select gpu, GPU type RTX 4090 (we need all the memory we can get and the 4090 allows for the highest requestable memory). I ran with a custom instance size of 5 cores with 16GB/core. You should not have to define a CUDA version.

In the notebook run the code below to make sure you are running on a gpu capable core. Sometimes gpu nodes have issues so be sure to check you are actually have acess to it.
```python
import torch
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print("Using:", device)
```

If you want to submit jobs to the cluster here is an example bash script. You should not have to define a specific node but it can be useful to do so if one of them is acting up.
    
    #!/bin/bash
    #SBATCH --job-name=gpu_CNN_LSTM
    #SBATCH --time=24:00:00
    #SBATCH --partition=gpu
    #SBATCH --gres=gpu:rtx4090:1
    #SBATCH --cpus-per-task=1
    #SBATCH --mem-per-cpu=80G
    #SBATCH --output=/storage/homefs/ts18c034/MasterKinase/models/LSTM/0.2dropout1layer0.0008lrcrossvalidation_fold4_shuffeled_output%j.txt 

    cd /storage/homefs/ts18c034/MasterKinase

    module load Anaconda3

    echo "Job started at $(date)"

    eval "$(conda shell.bash hook)"

    conda activate pytorchtest

    echo "starting python"

    exec python3 /storage/homefs/ts18c034/MasterKinase/scripts/some_script.py

The exec command from above makes the python script the head process, I found this to be useful when you have to deal with signals that get sent before your time runs out. Just be aware that if you are using the above, nothing below exec gets executed.
