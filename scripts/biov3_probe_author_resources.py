"""Probe official MAESTER resources for a verifiable publisher RDS link.

Failure/authorization challenge => NO-GO; no guessed link, label or download.
"""
from __future__ import annotations
import argparse, hashlib, json, re, urllib.parse, urllib.request, urllib.error
from html.parser import HTMLParser
from pathlib import Path
ROOT='https://vangalenlab.bwh.harvard.edu/resources/maester-2021/'
TARGET='TenX_CellLineMix_Seurat_Keep.rds'
class Links(HTMLParser):
    def __init__(self):
        super().__init__();self.hrefs=[]
    def handle_starttag(self,tag,attrs):
        if tag=='a':
            self.hrefs += [v for k,v in attrs if k=='href' and v]

def inspect_resource_index(raw:bytes)->dict:
    parser=Links();parser.feed(raw.decode('utf-8','replace'))
    found=[]
    for href in parser.hrefs:
        url=urllib.parse.urljoin(ROOT,href)
        parts=urllib.parse.urlparse(url)
        if parts.scheme=='https' and parts.hostname=='vangalenlab.bwh.harvard.edu' and parts.path.endswith('/'+TARGET):
            found.append(url)
    return {'root':ROOT,'index_sha256':hashlib.sha256(raw).hexdigest(),
            'verified_download_links':sorted(set(found)),
            'url_invented':False,
            'index_links_inspected':len(parser.hrefs)}

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path)
    args=p.parse_args()
    result={'source':'PUBLISHER_RESOURCE_INDEX','original_seurat_filename':TARGET,
            'original_rds_read':False,'barcode_to_line_verified':False,
            'k562_only_training_expression_written':False}
    try:
        request=urllib.request.Request(ROOT,headers={'User-Agent':'GCEV-BIOV3-author-identity-audit/0.1'})
        with urllib.request.urlopen(request,timeout=25) as response:
            raw=response.read(2_000_001)
        if len(raw)>2_000_000: raise ValueError('Publisher resource index too large')
        result.update(inspect_resource_index(raw))
        result['decision']='AUTHOR_RDS_LINK_DISCOVERED_BUT_NOT_YET_BARCODE_VALIDATED' if result['verified_download_links'] else 'AUTHOR_RDS_NOT_LINKED_IN_PUBLIC_INDEX'
    except (urllib.error.URLError, TimeoutError, ValueError) as exc:
        result['decision']='AUTHOR_RDS_NOT_PUBLICLY_RETRIEVABLE'
        result['access_error']=str(exc)[:400]
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
