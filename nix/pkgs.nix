{ postOverlays ? [] } :

let
  sources = import ./sources.nix;
  qchemOverlay = import sources.NixOS-QChem;
  nixpkgs = import sources.nixpkgs {
    overlays = [ qchemOverlay ] ++ postOverlays;
    allowUnfree = true;

    # See https://github.com/markuskowa/NixOS-QChem#configuration-via-nixpkgs
    qchem-config = {
      optAVX = true;
      useCuda = false;
    };
  };

in nixpkgs
