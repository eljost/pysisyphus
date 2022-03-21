let
  lock = builtins.fromJSON (builtins.readFile ../flake.lock);
  flakeCompat = import (fetchGit {
    url = "https://github.com/edolstra/flake-compat";
    rev = lock.nodes.flake-compat.locked.rev;
    ref = "master";
  }) { src = ./..; };

in flakeCompat.defaultNix
