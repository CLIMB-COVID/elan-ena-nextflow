# elan-ena

A Nextflow pipeline for smoothly transferring fresh consensus sequences to ENA via `webin-cli`.
This fork is maintained specifically for CLIMB-COVID assembly submissions to ENA.

## Parameters

### Required command line parameters

 | Name | Description |
 | ---- | ----------- |
 | `--study` | ENA study identifier (`PRJEB`) |
 | `--manifest` |  Assembly metadata manifest |
 | `--webin_jar` | Path to `webin-cli` JAR ([releases](https://github.com/enasequence/webin-cli/releases)) |
 | `--out` | Path to write successful accessions table |

### Required environment variables

Additionally, you will need to set the following parameters in your environment:

 | Name | Description |
 | ---- | ----------- |
 | `WEBIN_USER` | EMBL-EBI Webin username |
 | `WEBIN_PASS` | EMBL-EBI Webin password |

### Optional command line parameters

| Name | Description |
| ---- | ----------- |
| `--ascp` | Enable `ascp` transfer with `webin-cli` (`ascp` must be on your `PATH`) |
| `--test` | Enable `webin-cli` test mode |
| `--description` | Template for `DESCRIPTION` field (supports expansion of CSV `row` variables) |


## Invocation

```
export WEBIN_USER='Webin-00000'
export WEBIN_PASS='hunter2'
nextflow run climb-covid/elan-ena-nextflow -r stable \
    --study PRJEB00000 \
    --manifest /path/to/manifest.tsv \
    --webin_jar /path/to/webin-cli.jar \
    --out hoot.accessions.txt \
    --ascp --test
    --description 'SAMPLE:${-> row.central_sample_id}|RUN:${-> row.run_name}'
```

To update a local copy:

```
nextflow pull climb-covid/elan-ena-nextflow
```
