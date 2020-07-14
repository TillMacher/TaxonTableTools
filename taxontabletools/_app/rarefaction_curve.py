def rarefaction_curve_OTUs(TaXon_table_xlsx, replicates, error_style, rarefaction_ylim, path_to_outdirs):

    import random
    import PySimpleGUI as sg
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import plot, ion, show
    from pathlib import Path

    TaXon_table_file = Path(TaXon_table_xlsx)

    print("\n" + "Input file:", TaXon_table_file.name)

    output_file = Path(str(path_to_outdirs) + "/" + "Rarefaction_curves" + "/" + TaXon_table_file.name + "_rarefaction_otus.pdf")

    TaXon_table_xlsx = pd.ExcelFile(TaXon_table_xlsx)
    df = pd.read_excel(TaXon_table_xlsx, 'TaXon table', header=0)

    available_samples = df.columns.tolist()[10:]
    sample_dict_clean = {}

    # iterate through all available samples
    for sample in available_samples:
        # create a dict for the read numbers of the respective sample
        sample_dict = df[sample].to_dict()
        # create an empty list for keys to store
        key_list = []
        # iterate through the read table of the sample
        for key, value in sample_dict.items():
            # if the read number of the OTU is NOT ZERO
            if value != 0:
                # append the OTU number to the key list
                key_list.append(key)
        # add the key list of each sample to a clean dictionary
        sample_dict_clean[sample] = key_list

    # draw once for each sample
    number_of_draws = len(sample_dict_clean.keys())

    # dictionary to store the drawing results
    draw_dictionary = {}

    ############################################################################
    ## create the progress bar window
    layout = [[sg.Text('Progress bar')],
              [sg.ProgressBar(1000, orientation='h', size=(20, 20), key='progressbar')],
              [sg.Cancel()]]
    window_progress_bar = sg.Window('Progress bar', layout, keep_on_top=True)
    progress_bar = window_progress_bar['progressbar']
    progress_update = 0
    progress_increase = 1000 / replicates
    ############################################################################

    for n_reps in range(0, replicates):
        # store the original dictionary to start over again
        # a copy of the original dictionary is required, because the samples will be removed with each draw
        # thus for each replicate a new dictionary to draw from has to be created
        sample_dict_to_draw = dict(sample_dict_clean)

        OTU_list = []
        OTU_set = []

        for i in range(0, number_of_draws):
            # choose a random sample from the dictionary
            random_choice = random.choice(list(sample_dict_to_draw.keys()))
            # extract the OTU IDs from the chosen sample and add them to the already existing OTU IDs
            OTU_list = OTU_list + sample_dict_clean[random_choice]
            # create a unique set
            OTU_set = set(OTU_list)
            # number of OTUs
            n_OTUs = len(OTU_set)
            # now add the unique OTU list to the output dictionary
            # if the key is not in the dict, create a new entry (= OTU ID plus number of OTUs)
            if i not in draw_dictionary.keys():
                draw_dictionary[i] = [n_OTUs]
            # if the key already exists, calculate the sum of the already existing number of OTUs and the new number of OTUs
            else:
                # create a new list to store the current number of OTUs
                add_OTUs_list = draw_dictionary[i]
                add_OTUs_list.append(n_OTUs)
                draw_dictionary[i] = add_OTUs_list

            # remove the sample to draw only once
            sample_dict_to_draw.pop(random_choice)

        ############################################################################
        event, values = window_progress_bar.read(timeout=1)
        if event == 'Cancel'  or event is None:
            print('\nCancel')
            window_progress_bar.Close()
            raise RuntimeError
        # update bar with loop value +1 so that bar eventually reaches the maximum
        progress_update += progress_increase
        progress_bar.UpdateBar(progress_update)
        ############################################################################

    window_progress_bar.Close()

    # create a dict to store the average number of OTUs per draw
    rarefaction_dict_average, rarefaction_dict_stdef = {}, {}

    def average(lst):
        return sum(lst) / len(lst)

    # iterate through the draw_dictionary and calculate the average number of OTUs
    for key, value in draw_dictionary.items():
        average_OTUs = average(draw_dictionary[key])
        stdef_OTUs = np.std(draw_dictionary[key], dtype=np.float64)
        rarefaction_dict_average[key] = average_OTUs
        rarefaction_dict_stdef[key] = stdef_OTUs

    if error_style == "a":
        draws = [i+1 for i in rarefaction_dict_average.keys()] #list(rarefaction_dict_average.keys())
        n_OTUs = list(rarefaction_dict_average.values())
        error_bar = list(rarefaction_dict_stdef.values())
        plt.figure(figsize=(20, 10))
        if rarefaction_ylim != '':
            plt.ylim(0, int(rarefaction_ylim))
        plt.errorbar(draws, n_OTUs, error_bar, linewidth=0.8, color='blue', capsize=1.3, capthick=1, ecolor='lightgrey')
        plt.xticks(np.arange(1, len(draws)+1, step=1))
        plt.xlabel('# samples')
        plt.ylabel('# OTUs')
        plt.title('repetitions = ' + str(replicates))


    elif error_style == "b":
        draws = [i+1 for i in rarefaction_dict_average.keys()] #list(rarefaction_dict_average.keys())
        n_OTUs = list(rarefaction_dict_average.values())
        y = np.asarray(n_OTUs)
        error_bar = np.asarray(list(rarefaction_dict_stdef.values()))
        plt.figure(figsize=(20, 10))
        if rarefaction_ylim != '':
            plt.ylim(0, int(rarefaction_ylim))
        plt.plot(draws, n_OTUs)
        plt.fill_between(draws, y-error_bar, y+error_bar, alpha=0.1)
        plt.xticks(np.arange(1, len(draws)+1, step=1))
        plt.xlabel('# samples')
        plt.ylabel('# OTUs')
        plt.title('repetitions = ' + str(replicates))


    plt.show(block=False)
    answer = sg.PopupYesNo('Save figure?', keep_on_top=True)
    if answer == "Yes":
        plt.savefig(output_file, bbox_inches='tight')
        plt.close()
        closing_text = "\n" + "Rarefaction curves are found in: " + str(path_to_outdirs) + "/rarefaction_curves/"
        print(closing_text)
        sg.Popup(closing_text, title="Finished", keep_on_top=True)
    else:
        plt.close()


def rarefaction_curve_species(TaXon_table_xlsx, replicates, error_style, rarefaction_ylim, path_to_outdirs):

    import random
    import PySimpleGUI as sg
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import plot, ion, show
    from pathlib import Path

    TaXon_table_file = Path(TaXon_table_xlsx)

    print("\n" + "Input file:", TaXon_table_file.name)

    output_file = Path(str(path_to_outdirs) + "/" + "Rarefaction_curves" + "/" + TaXon_table_file.name + "_rarefaction_species.pdf")

    TaXon_table_xlsx = pd.ExcelFile(TaXon_table_xlsx)
    df = pd.read_excel(TaXon_table_xlsx, 'TaXon table', header=0)
    df = df.replace(np.nan,"nan")

    available_samples = df.columns.tolist()[10:]
    sample_dict_clean = {}

    # iterate through all available samples
    for sample in available_samples:
        # create a dict for the read numbers of the respective sample for each species
        sample_OTU_list = df[[sample, "Species"]].values.tolist()
        # select only the present Species
        sample_species_list = list(set([OTU[1] for OTU in sample_OTU_list if (OTU[0] != 0 and OTU[1] != "nan")]))
        # store the species in a dictionary
        sample_dict_clean[sample] = sample_species_list

    # draw once for each sample
    number_of_draws = len(sample_dict_clean.keys())

    # dictionary to store the drawing results
    draw_dictionary = {}

    ############################################################################
    ## create the progress bar window
    layout = [[sg.Text('Progress bar')],
              [sg.ProgressBar(1000, orientation='h', size=(20, 20), key='progressbar')],
              [sg.Cancel()]]
    window_progress_bar = sg.Window('Progress bar', layout, keep_on_top=True)
    progress_bar = window_progress_bar['progressbar']
    progress_update = 0
    progress_increase = 1000 / replicates
    ############################################################################

    for n_reps in range(0, replicates):
        # store the original dictionary to start over again
        # a copy of the original dictionary is required, because the samples will be removed with each draw
        # thus for each replicate a new dictionary to draw from has to be created
        sample_dict_to_draw = dict(sample_dict_clean)

        species_list = []
        species_set = []

        for i in range(0, number_of_draws):
            # choose a random sample from the dictionary
            random_choice = random.choice(list(sample_dict_to_draw.keys()))
            # extract the OTU IDs from the chosen sample and add them to the already existing OTU IDs
            species_list = species_list + sample_dict_clean[random_choice]
            # create a unique set
            species_set = set(species_list)
            # number of OTUs
            n_species = len(species_set)
            # now add the unique OTU list to the output dictionary
            # if the key is not in the dict, create a new entry (= OTU ID plus number of OTUs)
            if i not in draw_dictionary.keys():
                draw_dictionary[i] = [n_species]
            # if the key already exists, calculate the sum of the already existing number of OTUs and the new number of OTUs
            else:
                # create a new list to store the current number of OTUs
                add_species_list = draw_dictionary[i]
                add_species_list.append(n_species)
                draw_dictionary[i] = add_species_list

            # remove the sample to draw only once
            sample_dict_to_draw.pop(random_choice)

        ############################################################################
        event, values = window_progress_bar.read(timeout=1)
        if event == 'Cancel'  or event is None:
            print('\nCancel')
            window_progress_bar.Close()
            raise RuntimeError
        # update bar with loop value +1 so that bar eventually reaches the maximum
        progress_update += progress_increase
        progress_bar.UpdateBar(progress_update)
        ############################################################################

    window_progress_bar.Close()

    # create a dict to store the average number of OTUs per draw
    rarefaction_dict_average, rarefaction_dict_stdef = {}, {}

    def average(lst):
        return sum(lst) / len(lst)

    # iterate through the draw_dictionary and calculate the average number of OTUs
    for key, value in draw_dictionary.items():
        average_species = average(draw_dictionary[key])
        stdef_species = np.std(draw_dictionary[key], dtype=np.float64)
        rarefaction_dict_average[key] = average_species
        rarefaction_dict_stdef[key] = stdef_species

    if error_style == "a":
        draws = [i+1 for i in rarefaction_dict_average.keys()] #list(rarefaction_dict_average.keys())
        n_species = list(rarefaction_dict_average.values())
        error_bar = list(rarefaction_dict_stdef.values())
        plt.figure(figsize=(20, 10))
        if rarefaction_ylim != '':
            plt.ylim(0, int(rarefaction_ylim))
        plt.errorbar(draws, n_species, error_bar, linewidth=0.8, color='blue', capsize=1.3, capthick=1, ecolor='lightgrey')
        plt.xticks(np.arange(1, len(draws)+1, step=1))
        plt.xlabel('# samples')
        plt.ylabel('# species')
        plt.title('repetitions = ' + str(replicates))


    elif error_style == "b":
        draws = [i+1 for i in rarefaction_dict_average.keys()] #list(rarefaction_dict_average.keys())
        n_species = list(rarefaction_dict_average.values())
        y = np.asarray(n_species)
        error_bar = np.asarray(list(rarefaction_dict_stdef.values()))
        plt.figure(figsize=(20, 10))
        if rarefaction_ylim != '':
            plt.ylim(0, int(rarefaction_ylim))
        plt.plot(draws, n_species)
        plt.fill_between(draws, y-error_bar, y+error_bar, alpha=0.1)
        plt.xticks(np.arange(1, len(draws)+1, step=1))
        plt.xlabel('# samples')
        plt.ylabel('# species')
        plt.title('repetitions = ' + str(replicates))


    plt.show(block=False)
    answer = sg.PopupYesNo('Save figure?', keep_on_top=True)
    if answer == "Yes":
        plt.savefig(output_file, bbox_inches='tight')
        plt.close()
        closing_text = "\n" + "Rarefaction curves are found in: " + str(path_to_outdirs) + "/rarefaction_curves/"
        print(closing_text)
        sg.Popup(closing_text, title="Finished", keep_on_top=True)
    else:
        plt.close()
