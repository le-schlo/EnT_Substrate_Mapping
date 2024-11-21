import contextlib
import chemprop

home_dir = '/path/to/project_directory/'
ring_class = 'thiophene'
substitution = 'monosubstituted'

path_to_model = home_dir + 'Models/triplet_energy/triplet_energy.pt'
path_to_molecules = home_dir + f'Data/virtual_libraries/{ring_class}_{substitution}.csv'
path_for_predictions = home_dir + f'Data/et_predictions/predictions_{ring_class}_{substitution}.csv'

arguments = [
    '--test_path', path_to_molecules,
    '--preds_path', path_for_predictions,
    '--num_workers', '1',
    '--checkpoint_path', path_to_model,
    '--features_generator', 'rdkit_2d_normalized',
    '--no_features_scaling'
]

with contextlib.redirect_stdout(None):
    args = chemprop.args.PredictArgs().parse_args(arguments)
    preds = chemprop.train.make_predictions(args=args)
