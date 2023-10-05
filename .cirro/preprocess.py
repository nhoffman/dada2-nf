#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd


def setup_library(ds: PreprocessDataset) -> None:

    # Use the library param to set up the params appropriately
    alignment_options = {
        "SSU rRNA": {
            "library": "",
            "model": "data/SSU_rRNA_bacteria.cm",
            "strategy": "cmsearch"
        },
        "ssu-align-0.1.1-bacteria-0p1": {
            "library": "",
            "model": "data/ssu-align-0.1.1-bacteria-0p1.cm",
            "strategy": "cmsearch"
        },
        "Bacteria 16S Type Strains": {
            "library": "data/bacteria.16S.types.fasta.gz",
            "model": "",
            "strategy": "vsearch"
        },
        "ITS": {
            "library": "data/fungi.ITS.fna.gz",
            "model": "",
            "strategy": "vsearch"
        }
    }
    msg = f"Library option is not configured: {ds.params['library']}"
    assert alignment_options.get(ds.params['library']) is not None, msg
    ds.add_param(
        "alignment",
        alignment_options[ds.params['library']]
    )
    ds.remove_param('library')


def make_manifest(ds: PreprocessDataset) -> pd.DataFrame:

    # Make a manifest which will be passed to
    # the workflow as params.sample_information
    # Columns will be:
    #   sample, sampleid, datadir, R1, R2, I1, I2 (index reads optional)
    ds.logger.info(ds.files.to_csv(index=None))
    manifest = ds.files.assign(
        readType=lambda d: d['readType'].fillna("R"),
        readLabel=lambda d: d.apply(
            lambda r: f"{r['readType']}{int(r['read'])}",
            axis=1
        ),
        datadir=lambda d: d['file'].apply(
            lambda s: s.rsplit("/", 1)[0] + "/"
        ),
        file=lambda d: d['file'].apply(
            lambda s: s.rsplit("/", 1)[1]
        ),
        sampleid=lambda d: d['file'].apply(
            lambda fp: fp.split("_", 1)[0]
        )
    )

    ds.logger.info(manifest.to_csv(index=None))

    manifest = manifest.pivot(
        index=['sampleIndex', 'sample', 'sampleid', 'datadir'],
        columns='readLabel',
        values='file'
    )
    ds.logger.info(manifest.to_csv(index=None))

    manifest = manifest.reset_index(
    ).drop(
        columns=["sampleIndex"]
    ).rename(
        columns=dict(
            sample="sample_name"
        )
    )

    # The expected number of reads depends on the input_file_type param
    expected_reads = {
        'dual': ['R1', 'R2', 'I1', 'I2'],
        'single': ['R1', 'R2', 'I1'],
        'none': ['R1', 'R2'],
    }[ds.params.get("index_file_type")]

    first_columns = [
        "sampleid",
        "sample_name",
        "datadir"
    ] + expected_reads

    manifest = manifest.reindex(
        columns=first_columns + [
            cname
            for cname in manifest.columns.values
            if cname not in first_columns
        ]
    )

    ds.logger.info(manifest.to_csv(index=None))
    msg = "Error: All samples must be paired-end FASTQ"
    assert manifest.dropna().shape[0] == manifest.shape[0], msg

    # Add any other metadata uploaded by the user
    ds.logger.info("Adding user-provided metadata")
    ds.logger.info(ds.samplesheet.to_csv(index=None))
    manifest = manifest.assign(**{
        cname: manifest['sample_name'].apply(cvals.get)
        for cname, cvals in ds.samplesheet.set_index('sample').items()
    })
    ds.logger.info(manifest.to_csv(index=None))

    return manifest


def filter_manifest(
    ds: PreprocessDataset,
    manifest: pd.DataFrame
) -> pd.DataFrame:

    # Check and see if the user provided a manifest
    # If not, take no action
    manifest_fp = ds.params.get("manifest")
    if manifest_fp is None:
        return manifest

    # Read the user-provided manifest file
    ds.logger.info(f"Reading the user-provided manifest file ({manifest_fp})")
    if manifest_fp.endswith("csv"):
        ds.logger.info("Parsing CSV")
        user_manifest = pd.read_csv(manifest_fp)
    elif manifest_fp.endswith("xlsx"):
        ds.logger.info("Parsing XLSX")
        user_manifest = pd.read_excel(manifest_fp)
    else:
        raise Exception(f"Did not recognize file type: {manifest_fp}")
    ds.logger.info(user_manifest.to_csv(index=None))

    # Remove any whitespaces in the column headers
    user_manifest = user_manifest.rename(
        columns=lambda cn: cn.strip() if isinstance(cn, str) else cn
    )

    # Make sure that the user provided a 'sampleid' column
    msg = "User must include a 'sampleid' column in the manifest"
    assert 'sampleid' in user_manifest.columns.values, msg

    # Remove any whitespaces in the values of the manifest
    user_manifest = user_manifest.applymap(
        lambda v: v.strip() if isinstance(v, str) else v
    )

    ds.logger.info("Whitespace-trimmed manifest:")
    ds.logger.info(user_manifest.to_csv(index=None))

    # Make sure that each row in the user-provided manifest matches a sample
    name_map = dict()
    ds.logger.info("Checking for sample name mappings")
    for sampleid in user_manifest['sampleid'].values:

        # If the sampleid matches the sampleid column, take no action
        if sampleid in manifest["sampleid"].values:
            continue
        # If the sampleid matches a value in the sample_name column
        elif sampleid in manifest["sample_name"].values:

            # Get the sampleid value which matches
            corrected_sampleid = manifest.query(
                f"sample_name == '{sampleid}'"
            )[
                'sampleid'
            ].values[0]

            # Add it to the map
            name_map[sampleid] = corrected_sampleid

        # If there is no match
        else:
            # Raise an error
            raise Exception(f"Could not find any sample matching '{sampleid}'")

    # If there are any values which need to be corrected
    if len(name_map) > 0:

        # Rename the index values
        user_manifest = user_manifest.assign(
            sampleid=user_manifest['sampleid'].apply(
                lambda sn: name_map.get(sn, sn)
            )
        )

    # Now filter down the manifest to only include those samples
    manifest = manifest.loc[
        manifest['sampleid'].isin(user_manifest['sampleid'].values)
    ]
    ds.logger.info(f"Filtered the manifest down to {manifest.shape[0]:,} rows")

    # Add any other metadata provided by the user
    for cname, cvals in user_manifest.set_index('sampleid').items():
        ds.logger.info(f"Adding user-provided values for {cname}")
        manifest = manifest.assign(**{
            cname: manifest['sampleid'].apply(cvals.get)
        })

    ds.logger.info("Final filtered manifest")
    ds.logger.info(
        manifest.to_csv(index=None)
    )

    return manifest


if __name__ == "__main__":

    # Instantiate the Cirro dataset object
    ds = PreprocessDataset.from_running()

    # Set up the library information appropriately
    setup_library(ds)

    # Get the manifest with all files in the input dataset
    manifest = make_manifest(ds)

    # Filter to the list of samples specified by the user, if any
    manifest = filter_manifest(ds, manifest)

    # Write out the sample metadata
    manifest.to_csv("manifest.csv", index=None)
    ds.add_param("manifest", "manifest.csv", overwrite=True)

    # Set the fastq_list param to false
    ds.add_param("fastq_list", False)

    # Log the params
    ds.logger.info(ds.params)
