# -*- mode: python ; coding: utf-8 -*-


a = Analysis(
    ['scripts/GUI.py'],
    pathex=['.'],
    binaries=[],
    datas=[('scripts/File/icon10.png', '.'), ('scripts/File/meaw.ico', '.')],
    hiddenimports=['Bio.Align.fasta', 'logging', 'pkg_resources', 'setuptools'],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='PyCulturome',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon="scripts/File/meaw.ico"
)
