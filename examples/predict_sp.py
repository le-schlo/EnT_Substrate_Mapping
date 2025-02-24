import time
import torch
import pandas as pd

from EasyChemML.DataImport.DataImporter import DataImporter
from EasyChemML.DataImport.Module.CSV import CSV
from EasyChemML.Encoder import BertTokenizer
from EasyChemML.Environment import Environment
from EasyChemML.Encoder.impl_Tokenizer.SmilesTokenizer_SchwallerEtAll import SmilesTokenzier
from EasyChemML.Model.impl_Pytorch.Models.BERT.FP2MOL_BERT_Trans import FP2MOL_Bert

from rdkit import Chem

home_dir = '/path/to/project_directory/'
ring_class = 'thiophene'
substitution = 'monosubstituted'

file_path = home_dir + f'Data/virtual_libraries/{ring_class}_{substitution}.csv'
path_to_model = home_dir + 'Models/spin_population/spin_population.pt'
path_for_predictions = home_dir + f'Data/sp_predictions/predictions_{ring_class}_{substitution}.csv'

# ----------------------------------- Data Preprocessing -----------------------
device = "cuda:0"

###Do not change the following parameters unless you train a new model with different parameters
d_model = 512
heads = 4
N = 8
src_vocab_size =106
trg_vocab_size = 136
dropout = 0.1
max_seq_len = 100
###

model_object = FP2MOL_Bert(src_vocab_size, trg_vocab_size, N, heads, d_model, dropout, max_seq_len, device)
model_object.load_model_parameters(path_to_model)

env = Environment(WORKING_path_addRelativ='Output')

dataloader = {'data': CSV(file_path, columns=['smiles'])}
di = DataImporter(env)
bp = di.load_data_InNewBatchPartition(dataloader)

val_data_size = len(bp['data'])

print('Start BertTokenizer')
tokenizer = BertTokenizer()
tokenizer.convert(datatable=bp['data'], columns=['smiles'], n_jobs=4)

# ----------------------------------- Evaluation -----------------------

batch_size = 25
decoding_strategy = 'greedy'

s = SmilesTokenzier()
final_df=pd.DataFrame(data=None)

start_mol = 0
iter_steps = 0

start = time.time()
temp = start

iteration_per_chunk = int(val_data_size / batch_size)

for verbose, i in enumerate(range(int(iteration_per_chunk)+1)):

    end_mol = start_mol + batch_size

    true_smiSD = bp['data'][start_mol:end_mol]['smiles_ids']
    input_smi = bp['data'][start_mol:end_mol]['smiles_ids']
    input_smi = input_smi[:, 0:100]

    true_smi = bp['data'][start_mol:end_mol]['smiles']

    true_smiSD = torch.LongTensor(true_smiSD).to(torch.device(device))
    input_smi = torch.LongTensor(input_smi).to(torch.device(device))
    model_eval = model_object.fit_eval(input_smi, true_smiSD, method=decoding_strategy)
    outputs = model_eval.outputs
    start_mol = end_mol
    iter_steps += 1
    finished = False

    if decoding_strategy == 'beam_search':
        num_beams = 5
    else:
        num_beams = 1

    for example in range(batch_size):
        for beam in range(num_beams):
            counter = (example * num_beams) + beam
            pred_smi, pred_SD_string, DensityArray = s.getSmilesfromoutputwithSD(outputs[counter])
            
            df_example = pd.DataFrame({})
            df_example['Predicted_smiles'] = [pred_smi]
            df_example['Predicted_SD_string'] = [pred_SD_string]
            df_example['DensityArray'] = [DensityArray]
            df_example['True_smiles'] = [true_smi[example].decode("utf-8")]

            if pred_smi == str(true_smi[example].decode("utf-8")):
                df_example['Reconstrunction_error'] = [0]
            else:
                df_example['Reconstrunction_error'] = [1]

            final_df = pd.concat([final_df, df_example])
            if ((i*batch_size) + example+1) == val_data_size:
                finished = True
                break

        if finished:
            break
    
    print('----------------------------------')
    print(f'Done with {verbose+1} / {iteration_per_chunk} chunks')

    final_df.to_csv(path_for_predictions, index=False)
