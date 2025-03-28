import os
from helpers.j_helpers import run_cmd


"""
Module that generates the connectome and weighted matrices.
"""


def connectome_generation(paths, nthreads):
    """
    Generate the connectome matrix from the filtered tractogram and atlas parcels.
    Nodes and edges
    """

    tckgen_output = os.path.join(paths["tck_dir"], "tracks_10mio_minmax_restricted.tck")
    parcels_coreg = os.path.join(paths["atlas_dir"], "hcpmmp1_parcels_coreg.mif")
    connectome_invlength_csv = os.path.join(paths["atlas_dir"], "hcpmmp1_scale_invlength.csv")
    connectome_minmax_csv = os.path.join(paths["atlas_dir"], "hcpmmp1_minmax_rest.csv")
    assignments_minmax_csv = os.path.join(paths["atlas_dir"], "assignments_hcpmmp1_scale_length.csv")
    assignments_invlength_csv = os.path.join(paths["atlas_dir"], "assignments_hcpmmp1_scale_length.csv")
    sift2_output = os.path.join(paths["tck_dir"], "sift2weights.csv")


    """
    FF: As far as I understand, we dont use -scale_invnodevol nor -scale_invlenght,
    because f.e. -scale_invnodevol tells tck2connectome to scale each connection by 
    the inverse of the parcel volume, but we don't need that because we have the
    external streamline weights from SIFT2 now
    Are these combinable? 
    """
    if not os.path.exists(connectome_invlength_csv):
        run_cmd([
            "tck2connectome",
            "-tck_weights_in", sift2_output,
            tckgen_output,
            parcels_coreg,
            connectome_invlength_csv,
            "-out_assignment", assignments_invlength_csv,
            "-symmetric", "-zero_diagonal",
            "-scale_invlength",
            "-nthreads", str(nthreads),
            "-force"
        ])
    else:
        print(f"{connectome_invlength_csv} already exists")

    if not os.path.exists(connectome_invlength_csv):
        run_cmd([
            "tck2connectome",
            "-tck_weights_in", sift2_output,
            tckgen_output,
            parcels_coreg,
            connectome_minmax_csv,
            "-out_assignment", assignments_minmax_csv,
            "-symmetric", "-zero_diagonal",
            "-nthreads", str(nthreads),
            "-force"
        ])
    else:
        print(f"{connectome_minmax_csv} already exists")

    """
    mrview /home/nikita/Nikita_MRI/me/atlas/hcpmmp1_parcels_coreg.mif -connectome.init /home/nikita/Nikita_MRI/me/atlas/hcpmmp1_parcels_coreg.mif -connectome.load /home/nikita/Nikita_MRI/me/atlas/hcpmmp1.csv
    """

    # Representing nodes as cortical meshes
    mesh_obj = os.path.join(paths["atlas_dir"], "hcpmmp1_mesh.obj")
    run_cmd([
        "label2mesh", parcels_coreg, mesh_obj,
        "-force"
    ])

    # Constructing a representation of true edge routes with exemplars
    exemplars_invlength = os.path.join(paths["tck_dir"], "exemplars_invlength.tck")
    exemplars_minmax = os.path.join(paths["tck_dir"], "exemplars_minmax.tck")
    if not os.path.exists(exemplars_invlength):
        run_cmd([
            "connectome2tck", tckgen_output, assignments_invlength_csv,
            exemplars_invlength, "-tck_weights_in", sift2_output,
            "-exemplars", parcels_coreg,
            "-files", "single",
            "-nthreads", str(nthreads),
            "-force"
        ])
    if not os.path.exists(exemplars_minmax):
        run_cmd([
            "connectome2tck", tckgen_output, assignments_minmax_csv,
            exemplars_minmax, "-tck_weights_in", sift2_output,
            "-exemplars", parcels_coreg,
            "-files", "single",
            "-nthreads", str(nthreads),
            "-force"
        ])


def generate_weighted_connectome_matrices(paths, nthreads):
    """
    Generate dMRI metrics along streamlines and construct connectome matrices weighted by these metrics.

    This function performs the following steps:
      1. For each diffusion metric (FA, ADC, AD, RD), sample the underlying image along the streamlines
         to compute a mean value per streamline
      2. Generate connectome matrices weighted by each metric by scaling the streamline contributions
      3. Calculate the average streamline length between node pairs
    """

    tckgen_output = os.path.join(paths["tck_dir"], "tracks_10mio_minmax_restricted.tck")
    parcels_file = os.path.join(paths["atlas_dir"], "hcpmmp1_parcels_coreg.mif")
    weights_file = os.path.join(paths["tck_dir"], "sift2weights.csv")

    # Define output file paths for the streamline metric sampling.
    os.makedirs(paths["connectome_dir"], exist_ok=True)
    out_fa = os.path.join(paths["connectome_dir"], "tck_meanFA.csv")
    out_adc = os.path.join(paths["connectome_dir"], "tck_meanADC.csv")
    out_ad = os.path.join(paths["connectome_dir"], "tck_meanAD.csv")
    out_rd = os.path.join(paths["connectome_dir"], "tck_meanRD.csv")

    # Sample dMRI metrics along the streamlines (tcksample with "-stat_tck mean").
    run_cmd([
        "tcksample",
        tckgen_output,
        os.path.join(paths["dwi_dir"], "fa.mif"),
        out_fa,
        "-stat_tck", "mean",
        "-nthreads", str(nthreads),
        "-force"
    ])

    run_cmd([
        "tcksample",
        tckgen_output,
        os.path.join(paths["dwi_dir"], "adc.mif"),
        out_adc,
        "-stat_tck", "mean",
        "-nthreads", str(nthreads),
        "-force"
    ])

    run_cmd([
        "tcksample",
        tckgen_output,
        os.path.join(paths["dwi_dir"], "ad.mif"),
        out_ad,
        "-stat_tck", "mean",
        "-nthreads", str(nthreads),
        "-force"
    ])

    run_cmd([
        "tcksample",
        tckgen_output,
        os.path.join(paths["dwi_dir"], "rd.mif"),
        out_rd,
        "-stat_tck", "mean",
        "-nthreads", str(nthreads),
        "-force"
    ])

    # Generate connectome matrices weighted by each metric.
    # For FA-weighted connectome:
    run_cmd([
        "tck2connectome",
        tckgen_output,
        parcels_file,
        os.path.join(paths["connectome_dir"], "connectome_fa.csv"),
        "-tck_weights_in", weights_file,
        "-scale_file", out_fa,
        "-stat_edge", "mean",
        "-symmetric",
        "-zero_diagonal",
        "-nthreads", str(nthreads),
        "-force"
    ])

    # For ADC-weighted connectome:
    run_cmd([
        "tck2connectome",
        tckgen_output,
        parcels_file,
        os.path.join(paths["connectome_dir"], "connectome_adc.csv"),
        "-tck_weights_in", weights_file,
        "-scale_file", out_adc,
        "-stat_edge", "mean",
        "-symmetric",
        "-zero_diagonal",
        "-nthreads", str(nthreads),
        "-force"
    ])

    # For AD-weighted connectome:
    run_cmd([
        "tck2connectome",
        tckgen_output,
        parcels_file,
        os.path.join(paths["connectome_dir"], "connectome_ad.csv"),
        "-tck_weights_in", weights_file,
        "-scale_file", out_ad,
        "-stat_edge", "mean",
        "-symmetric",
        "-zero_diagonal",
        "-nthreads", str(nthreads),
        "-force"
    ])

    # For RD-weighted connectome:
    run_cmd([
        "tck2connectome",
        tckgen_output,
        parcels_file,
        os.path.join(paths["connectome_dir"], "connectome_rd.csv"),
        "-tck_weights_in", weights_file,
        "-scale_file", out_rd,
        "-stat_edge", "mean",
        "-symmetric",
        "-zero_diagonal",
        "-nthreads", str(nthreads),
        "-force"
    ])

    # Generate a connectome matrix of mean streamline lengths
    run_cmd([
        "tck2connectome",
        tckgen_output,
        parcels_file,
        os.path.join(paths["connectome_dir"], "meanlength.csv"),
        "-tck_weights_in", weights_file,
        "-scale_length",
        "-stat_edge", "mean",
        "-symmetric",
        "-zero_diagonal",
        "-nthreads", str(nthreads),
        "-force"
    ])

