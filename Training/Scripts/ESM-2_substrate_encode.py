import numpy as np
import pandas as pd
import torch
import seaborn as sns
import matplotlib.pyplot as plt
import torch.multiprocessing as mp
from transformers import AutoModel, AutoTokenizer
from sklearn.decomposition import PCA
from torch.utils.data import Dataset, DataLoader
from sklearn.preprocessing import LabelEncoder

import random
import math


model_src = "/storage/homefs/ts18c034/MasterKinase/ESM-2/esm2_t12_35M_UR50D"
tokenizer = AutoTokenizer.from_pretrained(model_src, from_tf=False)


class KinaseDataset(Dataset):
    def __init__(self, tokenized, labels, device):
        self.tokenized = tokenized
        self.labels = labels
        self.device = device

        # Convert labels to numerical format
        self.label_encoder = LabelEncoder()
        self.encoded_labels = torch.tensor(self.label_encoder.fit_transform(self.labels), dtype=torch.long)

    def __len__(self):
        return len(self.encoded_labels)

    def __getitem__(self, idx):
        tokenized_sample = {key: value[idx].to(self.device) for key, value in self.tokenized.items()}  
        label = self.encoded_labels[idx].to(self.device)

        return tokenized_sample, label

device = torch.device("cpu")

torch_model = AutoModel.from_pretrained(model_src, from_tf=False) #AutoModel already loads only the embedding layers and not the classifier head so we dont have to truncate the model.

torch_model.config.output_hidden_states = True


torch_model = torch_model.to(device) #"move" model to GPU if applicable


torch_model.eval()


all_embeddings = []
all_labels = []


df=pd.read_csv("/storage/homefs/ts18c034/MasterKinase/data/mouserat/deduplicated_mouse_fullSeq_centered.tsv",sep="\t")
unique_sequences=list(df["centered_sub_seq"].unique())

df["unique_sequence_mapping"]=None #contains information of which unique sequence it is: dont want to encode the same sequence over and over again
mapping = {seq: i for i, seq in enumerate(unique_sequences)}
df["unique_sequence_mapping"] = df["centered_sub_seq"].map(mapping)

#df.to_csv("/storage/homefs/ts18c034/MasterKinase/data/deduplicated_human_fullSeq_centered_filteredpos_uniqidx.tsv",index=False)
df.to_csv("/storage/homefs/ts18c034/MasterKinase/data/mouserat/deduplicated_mouse_fullSeq_centered_filteredpos_uniqidx.tsv",index=False)



tokenized = tokenizer(unique_sequences, return_tensors="pt",truncation=False) 
labels=np.arange(tokenized["input_ids"].shape[0])


dataset = KinaseDataset(tokenized, labels, device)
dataloader = DataLoader(dataset, batch_size=128, shuffle=False, num_workers=3)

with torch.no_grad():
    for batch in dataloader:
        tokenized_batch, labels_batch = batch  # Extract batch of tokenized sequences and labels
        
        # Pass batch through the model
        torch_output = torch_model(**tokenized_batch)

        #extract final layer
        final_layer_embeddings = torch_output.hidden_states[11]

        all_labels.append(labels_batch)
        all_embeddings.append(final_layer_embeddings)

for i, tensor in enumerate(all_embeddings):
    torch.save(tensor, f'/storage/homefs/ts18c034/MasterKinase/ESM-2/mouserat/substrateembedding/mouse_substrate_chunk{i}_35M.pt')

for j, label in enumerate(all_labels):
    torch.save(label, f'/storage/homefs/ts18c034/MasterKinase/ESM-2/mouserat/substrateembedding/mouse_substrate_label_chunk{j}_35M.pt')

