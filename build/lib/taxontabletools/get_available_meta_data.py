def get_available_meta_data(TaXon_table_xlsx, path_to_outdirs):

    import os
    import pandas as pd
    from pandas import DataFrame
    from pathlib import Path

    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    meta_data_xlsx = Path(str(path_to_outdirs) + "/" + "Meta_data_table" + "/" + TaXon_table_xlsx.stem + "_metadata.xlsx")

    if os.path.exists(meta_data_xlsx):

        meta_data_xlsx = pd.ExcelFile(meta_data_xlsx)
        data = pd.read_excel(meta_data_xlsx, header=0)

        available_meta_data = data.columns.tolist()[1:]
        available_meta_data = [str(i) for i in available_meta_data]

        return available_meta_data

    else:
        return False
