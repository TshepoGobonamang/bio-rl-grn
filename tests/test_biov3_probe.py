import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.biov3_probe_author_resources import inspect_resource_index

def test_publisher_source_link_must_be_verified_not_guessed():
    obj=inspect_resource_index(b'<a href="TenX_CellLineMix_Seurat_Keep.rds">source</a>')
    assert obj['verified_download_links']==['https://vangalenlab.bwh.harvard.edu/resources/maester-2021/TenX_CellLineMix_Seurat_Keep.rds']
    assert not obj['url_invented']

def test_external_hosts_and_similar_filenames_not_trusted():
    obj=inspect_resource_index(b'<a href="https://evil.example/TenX_CellLineMix_Seurat_Keep.rds">evil</a><a href="TenX_CellLineMix_Seurat_Keep.rds.fake">fake</a>')
    assert obj['verified_download_links']==[]
