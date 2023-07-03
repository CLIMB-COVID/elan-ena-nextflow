#!/usr/bin/env python
import sys
from bs4 import BeautifulSoup as bs
import re
import argparse
import datetime

from ocarina.client import Ocarina
from ocarina.util import get_config

parser = argparse.ArgumentParser()
parser.add_argument("--test", action="store_true", default=False)
parser.add_argument("webin_manifest")
parser.add_argument("webin_output_xml")
parser.add_argument("published_name")
args = parser.parse_args()

ocarina = Ocarina()
if args.test:
    ocarina.config = get_config(profile="test-service-outbound")
else:
    ocarina.config = get_config(profile="service-outbound")


def send_asm_accession(publish_group, accession_id, assemblyname, subm_date=None):
    if not accession_id:
        print("[FAIL] No accession provided for %s" % publish_group)
        return False

    if not subm_date:
        subm_date = datetime.datetime.now().strftime("%Y-%m-%d")

    success, obj = ocarina.api.put_accession(
        publish_group=publish_group,
        service="ENA-ASSEMBLY",
        accession=accession_id.strip(),
        accession2=assemblyname.strip(),
        public=True,
        public_date=subm_date,
    )
    if not success:
        print(
            "[FAIL] Failed to submit accession %s to PAG %s"
            % (accession_id, publish_group)
        )
    else:
        print("[OKAY] %s:%s" % (publish_group, accession_id))
    return success


# usage: webin_to_majora.py <webin_manifest> <webin_output_xml> <published_name>
published_name = args.published_name

assembly_name = None
for line in open(args.webin_manifest):
    k, v = line.strip().split(None, 1)
    if k == "ASSEMBLYNAME":
        assembly_name = v

if not assembly_name or not published_name:
    sys.exit(1)

fh = open(args.webin_output_xml)
soup = bs("".join(fh.readlines()), "xml")

try:
    erz = soup.findAll("ANALYSIS")[0]["accession"]
except:
    try:
        error_recipt = soup.findAll("ERROR")[0]
        erz = re.search("ERZ\d{5,9}", str(error_recipt))[0]
    except:
        sys.exit(2)

send_asm_accession(published_name, erz, assembly_name)

print("published_name", "assemblyname", "ena_assembly_id")
print(published_name, assembly_name, erz)
