import umap.umap_ as umap

def run_umap(data):
    reducer = umap.UMAP(
        n_components=3,
        n_neighbors=60,
        min_dist=0.15,
        metric='euclidean',
        random_state=0,
        n_jobs=1
    )
    return reducer.fit(data)
