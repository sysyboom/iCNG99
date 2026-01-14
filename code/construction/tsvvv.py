import pandas as pd


def extract_and_save(filename, output_file):
    data = pd.read_csv(filename, sep='\t')

    def extract_all_references(info, keyword):

        items = info.split(' ... ')
        matches = [item for item in items if keyword in item]
        if matches:
            return ' ... '.join(matches)  
        return "Not found"  

    data['bigg.reaction'] = data['xrefs'].apply(lambda x: extract_all_references(x, 'bigg.reaction'))
    data['seed.reaction'] = data['xrefs'].apply(lambda x: extract_all_references(x, 'seed.reaction'))
    data['metacyc.reaction'] = data['xrefs'].apply(lambda x: extract_all_references(x, 'metacyc.reaction'))
    #data['seed.compound'] = data['xrefs'].apply(lambda x: extract_all_references(x, 'seed.compound'))
    #data['reactome:'] = data['xrefs'].apply(lambda x: extract_all_references(x, 'reactome:'))

    output_data = data[['#query', 'mnx_id', 'bigg.reaction', 'seed.reaction','metacyc.reaction']]
    output_data.to_excel(output_file, index=False)

filename = '/mnt/NFS/fengch/new/annotationall/id-mapper43'  
output_file = '/mnt/NFS/fengch/new/annotationall/id-mapper43.xlsx'
extract_and_save(filename, output_file)
