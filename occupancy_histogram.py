import sqlite3
import os, sys
import pandas as pd
from iotbx.pdb import hierarchy
import libtbx.phil
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

# TODO Set up alternate simpler run script (not giant.jiffies.run_default)

from giant.jiffies import run_default

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
}
""", process_includes=True)


##############################################################
def append_multiple_databases():
    pass


def detect_repeat_soaks(params):
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
            cur.execute("SELECT CompoundSMILES, CrystalName, RefinementBoundConformation FROM mainTable WHERE "
                        "CompoundSMILES == ?", (compound[0],))
            repeat_xtals = cur.fetchall()

            # Output: CSV(s) that contains repeat soaks
            repeat_xtals_df = pd.DataFrame(repeat_xtals, columns=["CompoundSMILES", "CrystalName",
                                                                  "RefinementBoundConformation"])
            file_name = str(compound[1]) + '.csv'
            file_path = os.path.join(params.output.out_dir, file_name)
            repeat_xtals_df.to_csv(file_path)
            yield file_path, compound[1]

    # Close connection to the database
    cur.close()


def repeat_soak_has_bound_conformation(repeat_xtal_csv_path):
    """ Return Boolean based on whether repeat soak has any bound conformations"""
    repeat_xtal_df = pd.read_csv(repeat_xtal_csv_path, header=0)
    # Return False if repeat xtal dataframe is empty once NaN rows are removed
    return not (repeat_xtal_df.dropna().empty)


def get_pandda_lig_chain(pdb_path, params):
    # Open connection to sqlite database
    conn = sqlite3.connect(params.input.database_path)
    cur = conn.cursor()

    # Get crystal name given pdb path
    cur.execute("SELECT CrystalName FROM mainTable WHERE RefinementBoundConformation == ?", (pdb_path,))
    xtal_name = cur.fetchall()

    # Get ligand chain given crystal name
    cur.execute("SELECT PANDDA_site_ligand_chain FROM panddaTable WHERE CrystalName == ?", xtal_name[0])
    pandda_lig_chain = cur.fetchall()

    # Close connection to sqlite database
    cur.close()

    # convert from list of tuple, to value
    pandda_lig_chain = pandda_lig_chain[0]
    pandda_lig_chain = pandda_lig_chain[0]

    return pandda_lig_chain


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
        smiles = row['CompoundSMILES']
        occ_b_df = read_ligand_occupancy_b(pdb_path, params)

        yield occ_b_df


def all_lig_occ(repeat_xtal_csv_path, params):
    """ Get all ligands occupancy for a single repeat soak"""

    all_lig_occupancy = np.empty([1, 1])

    for occ_b_df in get_occupancy_b_df(repeat_xtal_csv_path, params):
        # Check whether the ligand has a single occupancy value
        # TODO Change to not rely on positional indexing of occupancy value, use column name
        if occ_b_df.apply(lambda x: x.nunique())[1] == 1:
            lig_occ = occ_b_df.loc("Occupancy")[0][1]
            print lig_occ
            all_lig_occupancy = np.append(all_lig_occ,lig_occ)
        else:
            print "Occupancy varies across ligand, histogram not currently generated"
            exit()

    return all_lig_occupancy


# TODO Add name of ligand to title of plot
def occupancy_histogram(compound, all_lig_occupancy, params):
    """Use Occupancy dataframes to generate histogram"""

    plt.hist(all_lig_occ, rwidth=0.75)
    plt.xlim(0, 1)
    plt.title("Occupancy Histogram: {}".format(compound))
    plt.xlabel("Occupancy")
    plt.ylabel("Frequency")
    plt.savefig("occ_hist.png", dpi=params.options.plot_dpi)
    plt.close()


def b_atom_plot(repeat_xtal_csv_path, params):
    """ Use Occupancy & B factor df to generate plots of B factor against atom name """

    for occ_b_df in get_occupancy_b_df(repeat_xtal_csv_path, params):
        cmap = matplotlib.cm.get_cmap('OrRd')
        plt.xlabel("Atoms")
        plt.ylabel("B Factor")
        ax = plt.gca()

        ax.xaxis.set_ticks(np.arange(len(occ_b_df['Atom'])))
        ax.xaxis.set_ticklabels(occ_b_df['Atom'], rotation=90)

        colour = cmap(occ_b_df["Occupancy"].mean())
        plt.plot(occ_b_df.B_factor, c=colour)
    plt.savefig("b_scatter.png", dpi=params.options.plot_dpi)
    plt.close()


def b_occ_scatter(repeat_xtal_csv_path, params):
    """ Plot Occupancy vs B factor for all ligands"""

    plt.xlabel("Occupancy")
    plt.ylabel("B Factor")

    all_occ_b_df = pd.DataFrame(columns=["Atom", "Occupancy", "B_factor"])

    for occ_b_df in get_occupancy_b_df(repeat_xtal_csv_path, params):
        plt.scatter(occ_b_df.Occupancy, occ_b_df.B_factor, marker="x")
        all_occ_b_df = all_occ_b_df.append(occ_b_df)

    # Linear fit         
    m, b = np.polyfit(all_occ_b_df.Occupancy, all_occ_b_df.B_factor, 1)
    plt.plot(all_occ_b_df.Occupancy, m * all_occ_b_df.Occupancy + b, '-')

    plt.savefig("b_occ.png", dpi=params.options.plot_dpi)


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


def html_report():
    """ Generate a HTML report based on repeat soak information"""
    pass


def run(params):
    # Make directory for output, if it doesn't exist
    if not os.path.exists(params.output.out_dir):
        os.mkdir(params.output.out_dir)

    # Loop over all compounds in database with more than params.options.min_repeat
    # If the compounds have bound conformations, generate plots to show occupancy
    for compound_csv, compound in detect_repeat_soaks(params):
        if repeat_soak_has_bound_conformation(compound_csv):
            all_ligand_occupancy = all_lig_occ(compound_csv, params)
            occupancy_histogram(compound, all_ligand_occupancy, params)
            b_atom_scatter_plot(compound_csv, params)
            b_occ_scatter(compound_csv, params)
        else:
            print "Repeat Soak with {} has no bound conformations".format(compound)


if __name__ == '__main__':
    run_default(
        run=run,
        master_phil=master_phil,
        args=sys.argv[1:],
        blank_arg_prepend=blank_arg_prepend,
        program=PROGRAM,
        description=DESCRIPTION
    )