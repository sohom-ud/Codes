'''
Loads plasma velocity moments from FPI.
'''

from pyspedas.projects.mms import fpi

def load_fpi_moms(trange, species="elc", probe="1", drate="brst", wipe=True):
    
    var_names = {
        "numberdensity": "N",
        "numberdensity_err": "N_err",
        "bulkv_gse": "v",
        "bulkv_err": "v_err",
        "bulkv_spintone_gse": "v_spin",
        "prestensor_gse": "Ptensor",
        "prestensor_err": "Ptensor_err" ,
        "temptensor_gse": "Temptensor"
    }

    data_dict = dict()

    dtype = "dis" if species == "ion" else "des"

    prefix = f"mms{probe}_{dtype}"

    varlist = dict()

    # Create list of variables
    for key, value in var_names.items():
        varlist[f'{prefix}_{key}_{drate}'] = f"{value}_{species}_{probe}"

    # Load FPI moment files
            
    try: 
        data = fpi(
            trange=trange,
            probe=probe,
            data_rate=drate,
            datatype=f"{dtype}-moms",
            center_measurement=True,
            varnames=list(varlist.keys()), 
            get_support_data=True,
            notplot=True
        )
    except:
        print("Moments data not found.")

    for key, value in data.items():

        skey = varlist[key]

        dim = len(data[key]['y'].shape) - 1

        data_dict[skey] = dict()

        data_dict[skey]['Epoch'] = data[key]['x']
        data_dict[skey]['Values'] = data[key]['y']
    
    if 'v' and 'v_spin' in var_names.values():

        data_dict[f'v_spincorr_{species}_{probe}'] = dict()
        data_dict[f'v_spincorr_{species}_{probe}']['Epoch'] = data_dict[f'v_{species}_{probe}']['Epoch']
        data_dict[f'v_spincorr_{species}_{probe}']['Values'] = data_dict[f'v_{species}_{probe}']['Values'] - data_dict[f'v_spin_{species}_{probe}']['Values']

    return data_dict
