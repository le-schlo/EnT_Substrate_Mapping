import base64
import os, dash, glob
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from umap import UMAP
from dash import dash_table
import plotly.express as px
from dash.dependencies import Input, Output
from dash import html
from dash import dcc

#--------------Define which data shall be visualised and set the corresponding paths-------------------#

ring_class = 'thiophene'
substitution = 'monosubstituted'
representation = 'ECFP' # Available options are ECFP, MACCS, MACCS+ECFP

home_dir = '/path/to/project_directory/'

picFolder_path = home_dir + f'Data/images/{ring_class}_{substitution}' #Important: The names of the files in the folder must be the index of the molecules in the data_path!
data_path = home_dir + f'Data/virtual_libraries/{ring_class}_{substitution}.csv'
e_t_path = home_dir + f'Data/et_predictions/predictions_{ring_class}_{substitution}.csv'

#--------------Load the data, compute representation and reduce dimensionality-------------------#

panda_dataframe = pd.read_csv(data_path)
panda_dataframe = pd.concat([panda_dataframe, pd.read_csv(e_t_path, usecols =['e_t'])], axis=1)
panda_dataframe['bins'] = pd.cut(panda_dataframe['e_t'], bins=np.arange(35,90,5), labels=np.arange(35,85,5))
panda_dataframe['IDS'] = panda_dataframe.index.to_list()

mol_list = [Chem.MolFromSmiles(smiles) for smiles in panda_dataframe['smiles']]

if representation == 'ECFP':
    fp_list = [AllChem.GetMorganFingerprintAsBitVect(molecule, radius=2, nBits=1024) for molecule in mol_list]
if representation == 'MACCS':
    fp_list = [MACCSkeys.GenMACCSKeys(molecule) for molecule in mol_list]
if representation == 'MACCS+ECFP':
    fp_list = [np.concatenate((MACCSkeys.GenMACCSKeys(molecule), AllChem.GetMorganFingerprintAsBitVect(molecule, radius=2, nBits=1024))) for molecule in mol_list]

FP_array = np.array(fp_list)
UMAPspace = UMAP(n_components=3,
                 n_neighbors=60,
                 min_dist=0.15,
                 metric='euclidean',
                 random_state=0,
                 n_jobs=1
                 ).fit(FP_array)

panda_dataframe['UMAP_1'] = UMAPspace.embedding_[:, 0]
panda_dataframe['UMAP_2'] = UMAPspace.embedding_[:, 1]
panda_dataframe['UMAP_3'] = UMAPspace.embedding_[:, 2]
data_array = panda_dataframe.to_numpy()

#--------------Load the images for easy analysis of the data. Differentiate three data situations-------------------#

globstatment = picFolder_path + '/*.jpeg'
imageFiles = glob.glob(globstatment)

if not imageFiles:                                                                                              # If no images are found in the folder, a dummy image is used
    print(f"No image files found in the folder: {picFolder_path} \n Using dummy image instead.")
    default_image = home_dir + r'\Data\images\default.png'
    pics = {}
    for id in panda_dataframe['IDS'].to_list():
        encoded_image = base64.b64encode(open(default_image, 'rb').read()).decode("utf-8")
        pics[id] = encoded_image

elif len(imageFiles) != len(data_array):                                                                        # If an inconsistend number of images is found (e.g, when the spin density map
    pics = {}                                                                                                   # is not available for all molecules), a dummy image is used in these cases
    print(f"Number of images in folder {picFolder_path} ({len(imageFiles)} images) does not match the number "
          f"of molecules in the dataset: {len(data_array)} \nUsing dummy image for missing images instead.")
    for id in panda_dataframe['IDS'].to_list():
        try:
            img_ = picFolder_path + f'/{id}.jpeg'
            encoded_image = base64.b64encode(open(img_, 'rb').read()).decode("utf-8")
        except:
            default_image = home_dir + r'\Data\images\default.png'
            encoded_image = base64.b64encode(open(default_image, 'rb').read()).decode("utf-8")
        pics[id] = encoded_image

else:                                                                                                            # If the number of images matches the number of molecules, the images are loaded normally
    pics = {}
    for i, file in enumerate(imageFiles):
        encoded_image = base64.b64encode(open(file, 'rb').read()).decode("utf-8")
        file_name = os.path.basename(file.split('.')[0])
        id = int(file_name)
        pics[id] = encoded_image

def generate_info_dataframe(mol_id):
    '''
    Function required for the interactive analysis. It generates a pandas dataframe containing
    the information about a molecule that is displayed in the table.
    '''
    index = int(np.where(data_transpose[:, 3] == mol_id)[0])
    mol_dataPoint = data_transpose[index]
    mol_SMILES = mol_dataPoint[0]
    mol_e_t = np.round(mol_dataPoint[1],2)

    dict = {}
    dict['Descriptors'] = ['ID', 'Smiles', 'Triplet energy [kcal/mol]']
    dict['Value'] = [mol_id, mol_SMILES, mol_e_t]
    return pd.DataFrame(dict)

#--------------Prepare the 3D plot-------------------#

data, data_transpose, mol_pics = data_array.transpose(), data_array, pics
plot_df = pd.DataFrame(data)
df = generate_info_dataframe(data_transpose[0][3])

# initiating the app
app = dash.Dash()  # defining the layout

fig = px.scatter_3d(x=data[4],
                    y=data[5],
                    z=data[6],
                    hover_name=data[3],
                    color=data[2],
                    color_discrete_map = {  35: '#e5f5e0',  # light green
                                            40: '#a1d99b',  # soft green
                                            45: '#74c476',  # medium green
                                            50: '#31a354',  # dark green
                                            55: '#006d2c',  # very dark green
                                            60: '#fcbba1',  # light red
                                            65: '#fc9272',  # medium red
                                            70: '#de2d26'},  # dark red
                    labels={'x': 'UMAP_1',
                            'y': 'UMAP_2',
                            'z': 'UMAP_3'},
                    height=900,
                    width=950)

fig.update_traces(marker=dict(size=4,
                              opacity=0.7,
                              line=dict(width=1,
                                        color='DarkSlateGrey')),
                  selector=dict(mode='markers'))

app.layout = html.Div([

    html.Div([
        dcc.Graph(id='DB-SMILES',
                  figure=fig
                  )],
        style={'width': '60%', 'display': 'inline-block', 'border-style': 'solid', 'vertical-align': 'middle'}),
    html.Div([
        html.Div([
            html.Img(id='mol_pic', src=f'data:image/png;base64,{mol_pics[0]}')
        ], id='mol_pic_div'),
        html.Div([
            dash_table.DataTable(
                id='mol_infotable',
                columns=[{"name": i, "id": i} for i in df.columns],
                data=df.to_dict('records')
            )], id='mol_infotable_div')], style={'width': '30%', 'display': 'inline-block',
                                                 'margin-left': '50px', 'border-style': 'solid',
                                                 'vertical-align': 'middle'})
])

@app.callback(Output('mol_infotable_div', 'children'), [Input('DB-SMILES', 'hoverData')])
def callback_stats(hoverData):
    if hoverData is None:
        mol_id = 0
    else:
        point = hoverData['points'][0]
        mol_id = int(point['hovertext'])

    out_infos = generate_info_dataframe(mol_id)

    dt = dash_table.DataTable(
        id='mol_infotable',
        columns=[{"name": i, "id": i} for i in out_infos.columns],
        data=out_infos.to_dict('records'),
        style_cell_conditional=[
            {'if': {'column_id': 'Descriptors'},
             'width': '30%'},
            {'if': {'column_id': 'Value'},
             'width': '70%'},
        ],
        style_data={
            'overflow': 'hidden',
            'textOverflow': 'ellipsis',
            'maxWidth': 0,
            'whiteSpace': 'normal',
        },
    )

    return dt

@app.callback(Output('mol_pic_div', 'children'), [Input('DB-SMILES', 'hoverData')])
def callback_stats(hoverData):
    if hoverData is None:
        mol_id = int(0)
        return html.Img(
            id='mol_pic',
            src=f'data:image/png;base64,{mol_pics[mol_id]}',
            style={'width': '85%', 'height': 'auto', 'display': 'block', 'margin': 'auto'}
        )
    else:
        point = hoverData['points'][0]
        mol_id = int(point['hovertext'])
        return html.Img(
            id='mol_pic',
            src=f'data:image/png;base64,{mol_pics[mol_id]}',
            style={'width': '85%', 'height': 'auto', 'display': 'block', 'margin': 'auto'}
        )

if __name__ == '__main__':
    app.run_server(host='0.0.0.0')
