import sqlite3
import os
import sys
import pandas as pd
from iotbx.pdb import hierarchy
import libtbx.phil
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

# TODO Set up alternate simpler run script (not giant.jiffies.run_default)

from giant.jiffies import run_default
import giant.xray.edstats as ed
##############################################################

PROGRAM = 'occupancy_report'
DESCRIPTION = """
    Read in XChem explorer database(s), and file structure. 
    Extract repeat soaking examples, produce a report of interesting plots.
"""
blank_arg_prepend = {}

##############################################################
master_phil = libtbx.phil.parse("""
input{
    database_path = None
        .type = str
    alt_database_path = None
        .type = str
}
output{
    out_dir = "./Output"
        .type = str
}
options{
    min_repeats = 5
        .type = int 
    plot_dpi = 300
        .type = float
    resolution_limit = 4.0
        .type = float
    overwrite = False
        .type = bool
}
""", process_includes=True)


##############################################################

# TODO Write a way to copy files out from database/ file structure to be minimal testing set on Zambezi & home PC

def append_multiple_databases():
    pass


# If we use Refinement bound conformation, then we are ignoring the pandda refinement into a merged model. This is
# ok for checking occupancy, but problematic as soon as we want to run things like edstats on the model. If we use
# refinement pdb then we need to deal with the multiple copies of the ligand problem.

def detect_repeat_soaks(protein_name, params):
    """ From XCE SQL database detect which ligands are repeated. Output to CSV."""

    # Open connection to sqlite database
    conn = sqlite3.connect(params.input.database_path)
    cur = conn.cursor()

    # Get SMILES with at least the minimum number of repeats
    min_repeat_tuple = (params.options.min_repeats,)
    cur.execute('SELECT CompoundSMILES,CompoundCode, COUNT(CompoundSMILES) FROM '
                'mainTable GROUP BY CompoundSMILES HAVING (COUNT(CompoundSMILES) >=?)', min_repeat_tuple)
    repeat_smiles = cur.fetchall()
    # Read mainTable
    # Input: Database path via params
    # Options: Flag for minimum number of repeats.
    # Alter to tuple for input into C execute

    # Create a CSV for each repeat
    for compound in repeat_smiles:
        if 'None' == compound[0]:
            continue
        elif compound[0] is u'':
            continue
        else:
            cur.execute("SELECT CompoundSMILES, CrystalName, RefinementBoundConformation, RefinementPDB_latest,"
                        "RefinementMTZ_latest  FROM mainTable WHERE "
                        "CompoundSMILES == ?", (compound[0],))
            repeat_xtals = cur.fetchall()

            # Output: CSV(s) that contains repeat soaks
            repeat_xtals_df = pd.DataFrame(repeat_xtals, columns=["CompoundSMILES", "CrystalName",
                                                                  "RefinementBoundConformation",
                                                                  "RefinementPDB_latest",
                                                                  "RefinementMTZ_latest"])
            file_name = str(compound[1]) + '.csv'
            file_path = os.path.join(params.output.out_dir, protein_name,  file_name)
            repeat_xtals_df.to_csv(file_path)

            yield file_path, compound[1]

    # Close connection to the database
    cur.close()


def get_protein_name(params):
    """ Get protein name from database. Check database has only one protein in it"""

    # Open connection to sqlite database
    conn = sqlite3.connect(params.input.database_path)
    cur = conn.cursor()

    # Get Protein Name
    cur.execute("SELECT DISTINCT ProteinName FROM mainTable")
    clean_protein_list = []
    protein_names = cur.fetchall()
    for protein in protein_names:
        if protein[0] is None:
            continue
        elif protein[0] is u'':
            continue
        elif protein[0] in u'None':
            continue
        else:
            clean_protein_list.append(protein[0])
    if len(clean_protein_list) == 1:
        return clean_protein_list[0]
    else:
        print "Mulitple proteins in single datasource, Quitting as this needs handling seperately"
        # TODO Change to better error handling, and skip this database
        os._exit(1)

    cur.close()


def repeat_soak_has_bound_conformation(repeat_xtal_csv_path):
    """ Return Boolean based on whether repeat soak has any bound conformations"""

    repeat_xtal_df = pd.read_csv(repeat_xtal_csv_path, header=0)

    # Check that pdb files exist
    for pdb in repeat_xtal_df.dropna()["RefinementBoundConformation"]:
        os.path.exists(pdb)
        print "PDB file at {} does not exist".format(pdb)
        return False

    for pdb in repeat_xtal_df.dropna()["RefinementPDB_latest"]:
        os.path.exists(pdb)
        print "PDB file at {} does not exist".format(pdb)
        return False

    # Return False if repeat xtal dataframe is empty once NaN rows are removed
    return not repeat_xtal_df.dropna().empty


def get_pandda_lig_chain(pdb_path, params):
    # Open connection to sqlite database
    conn = sqlite3.connect(params.input.database_path)
    cur = conn.cursor()

    # Get crystal name given pdb path (can use split.bound or refine)
    cur.execute("SELECT CrystalName FROM mainTable WHERE RefinementBoundConformation == ?", (pdb_path,))
    xtal_name = cur.fetchall()
    if not xtal_name:
        cur.execute("SELECT CrystalName FROM mainTable WHERE RefinementPDB_latest == ?", (pdb_path,))
        xtal_name = cur.fetchall()

    # Get ligand chain given crystal name
    cur.execute("SELECT PANDDA_site_ligand_chain FROM panddaTable WHERE CrystalName == ?", xtal_name[0])
    pandda_lig_chain = cur.fetchall()
    # Close connection to sqlite database
    cur.close()

    if pandda_lig_chain:
        # convert from list of tuple, to value
        pandda_lig_chain = pandda_lig_chain[0]
        pandda_lig_chain = pandda_lig_chain[0]
        return pandda_lig_chain
    else:
        print "{} at {} does not appear to have a pandda ligand chain,".format(xtal_name, pdb_path)
        print "perhaps your datasource is pointing to the wrong file?"


def read_ligand_occupancy_b(pdb_path, params):
    """Extract occupancy and B factor of ligand of interest from one PDB file into a dataframe"""

    # Input: A PDB structure. XCE database via params
    # Options: Read the surrounding structure as well as the ligand. Angstrom distance?
    # Output: Occupancy for ligand in supplied pdb. Dict including chain & altloc?

    # Get ligand chain that is associated with Event
    pandda_lig_chain = get_pandda_lig_chain(pdb_path, params)

    # Read in single PDB file
    pdb_in = hierarchy.input(file_name=pdb_path)
    sel_cache = pdb_in.hierarchy.atom_selection_cache()
    lig_sel = sel_cache.selection("chain {}".format(pandda_lig_chain))
    lig_hierarchy = pdb_in.hierarchy.select(lig_sel)

    lig_occ_b = []
    # Get occupancy & B factor of ligand
    for model in lig_hierarchy.models():
        for chain in model.chains():
            for rg in chain.residue_groups():
                for ag in rg.atom_groups():
                    for atom in ag.atoms():
                        lig_occ_b.append([atom.name, atom.occ, atom.b])
    occ_b_df = pd.DataFrame(lig_occ_b, columns=["Atom", "Occupancy", "B_factor"])

    return occ_b_df


def get_occupancy_b_df(repeat_xtal_csv_path, params):
    """ Generate dataframe for each pdb that contains the B factor and occupancy of the ligand"""

    repeat_xtal_df = pd.read_csv(repeat_xtal_csv_path, header=0)
    repeat_xtal_clean = repeat_xtal_df.dropna()
    for index, row in repeat_xtal_clean.iterrows():
        pdb_path = row['RefinementBoundConformation']
        occ_b_df = read_ligand_occupancy_b(pdb_path, params)

        yield occ_b_df


def get_pdb_mtz_name(repeat_xtal_csv_path, params):
    """ Get refienement pdb path and mtz path from csv"""

    repeat_xtal_df = pd.read_csv(repeat_xtal_csv_path, header=0)
    repeat_xtal_clean = repeat_xtal_df.dropna()
    for index, row in repeat_xtal_clean.iterrows():
        pdb_path = row['RefinementPDB_latest']
        mtz_path = row['RefinementMTZ_latest']
        xtal_name = row['CrystalName']

        yield pdb_path,mtz_path, xtal_name


def all_lig_occ(repeat_xtal_csv_path, params):
    """ Get all ligands occupancy for a single repeat soak"""

    all_lig_occupancy = []

    for occ_b_df in get_occupancy_b_df(repeat_xtal_csv_path, params):
        # Check whether the ligand has a single occupancy value
        # TODO Change to not rely on positional indexing of occupancy value, use column name
        if occ_b_df.apply(lambda x: x.nunique())[1] == 1:
            lig_occ = occ_b_df.loc("Occupancy")[0][1]
            all_lig_occupancy.append(lig_occ)
        else:
            print "Occupancy varies across ligand, histogram not currently generated"
            exit()

    return all_lig_occupancy


def occupancy_histogram(protein_name, compound, all_lig_occupancy, params):
    """Use Occupancy dataframes to generate histogram"""

    if len(all_lig_occupancy) == 1:
        print "Only one soak, no histogram generated"
        return

    plt.hist(all_lig_occupancy, rwidth=0.75)
    plt.xlim(0, 1)
    plt.title("Occupancy Histogram: {} : {}".format(protein_name, compound))
    plt.xlabel("Occupancy")
    plt.ylabel("Frequency")

    file_path = os.path.join(params.output.out_dir, protein_name, compound, "occ_hist.png")

    plt.savefig(file_path, dpi=params.options.plot_dpi)
    plt.close()


# TODO Add Colourbar
def b_atom_plot(protein_name, compound, repeat_xtal_csv_path, params):
    """ Use Occupancy & B factor df to generate plots of B factor against atom name """

    fig, ax = plt.subplots()
    plt.title("{} {}".format(protein_name, compound))
    plt.xlabel("Atoms")
    plt.ylabel("B Factor")
    cmap = matplotlib.cm.get_cmap('OrRd')

    for occ_b_df in get_occupancy_b_df(repeat_xtal_csv_path, params):
        ax.xaxis.set_ticks(np.arange(len(occ_b_df['Atom'])))
        ax.xaxis.set_ticklabels(occ_b_df['Atom'], rotation=90)

        colour = cmap(occ_b_df["Occupancy"].mean())
        plt.plot(occ_b_df.B_factor, c=colour)

    file_path = os.path.join(params.output.out_dir, protein_name, compound, "b_scatter.png")

    plt.savefig(file_path, dpi=params.options.plot_dpi)
    plt.close()


def b_occ_scatter(protein_name, compound, repeat_xtal_csv_path, params):
    """ Plot Occupancy vs B factor for all ligands"""

    plt.title("{} {}".format(protein_name, compound))
    plt.xlabel("Occupancy")
    plt.ylabel("B Factor")
    plt.xlim(0, 1)

    all_occ_b_df = pd.DataFrame(columns=["Atom", "Occupancy", "B_factor"])

    for occ_b_df in get_occupancy_b_df(repeat_xtal_csv_path, params):
        plt.scatter(occ_b_df.Occupancy, occ_b_df.B_factor, marker="x")
        all_occ_b_df = all_occ_b_df.append(occ_b_df)

    # Linear fit         
    m, b = np.polyfit(all_occ_b_df.Occupancy, all_occ_b_df.B_factor, 1)
    plt.plot(all_occ_b_df.Occupancy, m * all_occ_b_df.Occupancy + b, '-')

    file_path = os.path.join(params.output.out_dir, protein_name, compound, "b_occ.png")

    plt.savefig(file_path, dpi=params.options.plot_dpi)
    plt.close()


# TODO Understand how occupancy convergence is measured in refmac refinements.
# TODO Dependency to have used refmac?
def occupancy_convergence_plots():
    pass


# TODO Decide whether PLIF style plot is needed
def occupancy_plif_plot():
    """ Produce a PLIF plot that also displays occuapncy & """
    pass


def occupancy_soak_time():
    """ Generate a plot that compares occupancy to soak time for all hits """
    pass


# TODO Write Exhaustive search code
# TODO Decide how to compare exhaustive search and refinement
def refinement_vs_exhaustive():
    pass


# TODO Plots comparing the attempted repeat soak with reason for failure
def get_soaked_crystals(compound,params):

    # Open connection to sqlite database
    conn = sqlite3.connect(params.input.database_path)
    cur = conn.cursor()

    # Get number of soak attempts
    cur.execute("SELECT crystalName From mainTable WHERE compoundCode == ?", (compound,))
    soak_attempts = cur.fetchall()
    num_soak_attempts = len(soak_attempts)
    num_failed_soaks = 0
    for crystal in soak_attempts:
        if crystal[0] == u'':
            num_failed_soaks += 1
    num_soaks = num_soak_attempts - num_failed_soaks

    # Get number of diffracting xtals
    cur.execute("SELECT DimpleResolutionHigh FROM mainTable WHERE compoundCode == ? AND" 
                " DimpleResolutionHigh IS NOT NULL AND DimpleResolutionHigh < ?",
                (compound,params.options.resolution_limit,))
    diffracting_xtals = cur.fetchall()
    num_diffracting_xtals = len(diffracting_xtals)

    # Get number of bound conformations & events
    cur.execute("SELECT mainTable.RefinementBoundConformation FROM panddaTable, mainTable "
                 "WHERE mainTable.CrystalName = panddaTable.CrystalName AND compoundCode == ?", (compound,))
    events = cur.fetchall()
    num_events = len(events)
    num_failed_events = 0
    if num_events != 0:
        for event in events:
            if event[0] == None:
                num_failed_events +=1
    num_bound_events = num_events-num_failed_events

    # Close connection to the database
    cur.close()

    return [num_soak_attempts, num_soaks, num_diffracting_xtals, num_events, num_bound_events]


def soak_failure_bar_chart(protein_name,compound,params):

    repeat_soak_numbers = get_soaked_crystals(compound,params)
    N = len(repeat_soak_numbers)
    x = range(N)
    width = 0.75
    plt.bar(x,repeat_soak_numbers,width, color ="green")

    labels = ("Soak Attempts", "Successful Soaks", "Xtals diffracting to < {} Angstrom".format(compound),
              "Xtals with events", "Xtals with bound conformations")

    fig, ax = plt.subplots()

    plt.title("{} {}".format(protein_name, compound))
    plt.ylabel("Number of Xtals")
    ax.set_xticklabels(labels)

    file_path = os.path.join(params.output.out_dir, protein_name, compound, "soak_failure.png")
    plt.savefig(file_path, dpi=params.options.plot_dpi)
    plt.close()


# TODO Decide whether report is needed
def html_report():
    """ Generate a HTML report based on repeat soak information"""
    pass


def edstats(mtz_file_path, pdb_file_path, edstats_csv_path):
    """ Run Edstats on a single file"""

    # Running Edstats
    edstats, summary = ed.score_file_with_edstats(mtz_file_path,pdb_file_path)
    edstats.scores.to_csv(edstats_csv_path)

    yield edstats.scores

def residue_plot_edstats(edstats_scores, key, ylabel, protein_name, compound, xtal_name, params):

    # TODO Replace with more general form, splitting out plotting to a seperate fucntion
    # TODO Replace "A" with pulling out chains of the protein
    # Splitting RSR score into required chain
    main_chain = edstats_scores.loc[key, (slice(None), "A", slice(None), slice(None))]
    ordered_chain = main_chain.sort_index(level=2)

    # Plotting of edstats row vs residue for one chain
    fig = plt.figure()
    plt.plot(ordered_chain.index.get_level_values(2).values, ordered_chain.values)
    plt.ylabel(ylabel)
    plt.xlabel("Residue Number")
    img_file_path = os.path.join(params.output.out_dir, protein_name, compound,
                                 "{}_{}".format(ylabel.replace(' ','_'),xtal_name))
    plt.savefig(img_file_path, dpi = params.options.plot_dpi)
    plt.close()


def collate_edstats(repeat_xtal_csv_path, compound_name, protein_name, params):
    """ Get edstats results from all xtals in a repeat soak"""

    all_edstats_scores_list = []
    xtal_names = []
    lig_check = set()

    for pdb, mtz, xtal_name in get_pdb_mtz_name(repeat_xtal_csv_path, params):
        edstats_csv_path = os.path.join(params.output.out_dir, protein_name,
                                        compound_name, "edstat_{}.csv".format(xtal_name))
        for edstats_scores in edstats(mtz, pdb, edstats_csv_path):
            all_edstats_scores_list.append(edstats_scores)

            # TODO Move these plots to a loop over the ouputted csv file
            residue_plot_edstats(edstats_scores, 'Ra', 'RSR',protein_name, compound_name, xtal_name, params)
            residue_plot_edstats(edstats_scores, 'BAa', 'Mean residue B factor', protein_name, compound_name, xtal_name,params)
            residue_plot_edstats(edstats_scores, 'CCPa', 'RSCC (population)',protein_name, compound_name, xtal_name, params)
            residue_plot_edstats(edstats_scores, 'CCSa', 'RSCC (sample)',protein_name, compound_name, xtal_name, params)
            residue_plot_edstats(edstats_scores, 'ZDa', 'RSZD',protein_name, compound_name, xtal_name, params)
            residue_plot_edstats(edstats_scores, 'ZOa', 'RSZO',protein_name, compound_name, xtal_name, params)
            residue_plot_edstats(edstats_scores, 'ZD+a', 'RSZD+',protein_name, compound_name, xtal_name, params)
            residue_plot_edstats(edstats_scores, 'ZD+a', 'RSZD-',protein_name, compound_name, xtal_name, params)

        xtal_names.append(xtal_name)
        pandda_lig_chain = get_pandda_lig_chain(pdb,params)
        lig_check.add(pandda_lig_chain)

    if len(lig_check) == 1:
        pass
    elif len(lig_check) == 0:
        print " No ligand chains are present?"
        return
    else:
        print lig_check
        print "Ligand is not in same chain for all pdbs"
        return

    try:
        all_edstats_scores = pd.concat(all_edstats_scores_list, axis = 1, keys = xtal_names)
        file_path = os.path.join(params.output.out_dir, protein_name, compound_name,
                                 "all_edstats_{}.csv".format(compound_name))
        all_edstats_scores.to_csv(file_path)
        return lig_check, all_edstats_scores

    # TODO Make this handle the correct exception, rather than generic
    except Exception:
        print "No edstats calculated when bound conformation does not exist"


# TODO make general to plot any parameter sorted in edstats
def plot_mean_RSR(all_edstats_scores, protein_name, compound, pandda_lig_chain, params):
    """ Given a dataframe containing multiple xtal edstats scores return a mean RSR plot"""

    # TODO Replace selection of "A" chain by something that find protein chains.

    img_file_path = os.path.join(params.output.out_dir, protein_name, compound, "RSR_Mean")

    if all_edstats_scores is not None:

        RSR_chain = all_edstats_scores.loc['Ra', (slice(None),slice(None), "A", slice(None), slice(None))]
        RSR_chain = RSR_chain.sort_index(level=2)

        fig = plt.figure()
        plt.plot(np.sort(RSR_chain.index.get_level_values(level = 3).unique()), RSR_chain.groupby(level = 3).mean().values)
        plt.fill_between(np.sort(RSR_chain.index.get_level_values(level = 3).unique()),
                         RSR_chain.groupby(level=3).mean().values - RSR_chain.groupby(level = 3).std().values,
                         RSR_chain.groupby(level=3).mean().values + RSR_chain.groupby(level = 3).std().values,
                         alpha=0.2)
        plt.ylabel("RSR")
        plt.xlabel("Residue Number")



        lig_chain = str(pandda_lig_chain.pop())

        # If ligand chain cannot be pulled out from the edstats scores, then finish the plot without the ligand
        try:
            RSR_LIG = all_edstats_scores.loc['Ra', (slice(None),slice(None), lig_chain, slice(None), slice(None))]
        except Exception as error_inst:
            print "Ligand extaction failing: {}".format(error_inst)
            plt.savefig(img_file_path, dpi=params.options.plot_dpi)
            plt.close()
            return

        # TODO find way to plot near the closest residue.
        # Also consider a plot which takes the residues within nAng of ligand?
        plt.errorbar(x = 250,y = RSR_LIG.mean(), yerr = RSR_LIG.std(), fmt='o' )
        plt.savefig(img_file_path, dpi = params.options.plot_dpi)
        plt.close()


def run(params):
    # Make directory for output, if it doesn't exist
    if not os.path.exists(params.output.out_dir):
        os.mkdir(params.output.out_dir)

    protein_name = get_protein_name(params)

    # Make directory for each protein
    if not os.path.exists(os.path.join(params.output.out_dir, protein_name)):
        os.mkdir(os.path.join(params.output.out_dir, protein_name))

    # Loop over all compounds in database with more than params.options.min_repeat
    # If the compounds have bound conformations, generate plots to show occupancy

    # TODO Add ability to loop over multiple database files
    # TODO Add flag to not rerun data already collected
    # TODO Add remote connection to diamond, so can be run/ developed without diamond access
    # TODO Add ability to detect whether file exist, and and a overwrite parameter to make stuff anyways

    for compound_csv, compound in detect_repeat_soaks(protein_name, params):

        # TODO fix If/else wrapper to loading csv, the multilevel index is not loading well from the csv
        edstats_path = os.path.join(params.output.out_dir, protein_name, compound,"all_edstats_{}.csv".format(compound))

        try:
            pandda_lig_chain, all_edstats_scores = collate_edstats(compound_csv, compound, protein_name, params)
        except TypeError as error_inst:
            print "A type error has occured \n {} , ".format(error_inst)
            print "likely due to collate_edstats failing, perhaps due to ligand chain not being supplied"
            continue

        plot_mean_RSR(all_edstats_scores, protein_name, compound, pandda_lig_chain, params)

        """"
        print all_edstats_scores
        if os.path.exists(edstats_path):
            all_edstats_scores = pd.read_csv(edstats_path, index_col=[0], header=[0,1,2,3,4])
        else:
            all_edstats_scores = collate_edstats(compound_csv, compound, protein_name, params)
        """

        """
        try:
            plot_mean_RSR(all_edstats_scores)
        except:
            continue
        """
        """"# Make directory for each compound
        if not os.path.exists(os.path.join(params.output.out_dir, protein_name, compound)):
            os.mkdir(os.path.join(params.output.out_dir, protein_name, compound))

        if repeat_soak_has_bound_conformation(compound_csv):
            print "Generating plots for {} with {} bound".format(protein_name, compound)
            soak_failure_bar_chart(protein_name, compound, params)
            all_ligand_occupancy = all_lig_occ(compound_csv, params)
            occupancy_histogram(protein_name, compound, all_ligand_occupancy, params)
            b_atom_plot(protein_name, compound, compound_csv, params)
            b_occ_scatter(protein_name, compound, compound_csv, params)
        else:
            print "{} Repeat Soak with {} has no bound conformations".format(protein_name, compound)
        """
if __name__ == '__main__':
    run_default(
        run=run,
        master_phil=master_phil,
        args=sys.argv[1:],
        blank_arg_prepend=blank_arg_prepend,
        program=PROGRAM,
        description=DESCRIPTION
    )
