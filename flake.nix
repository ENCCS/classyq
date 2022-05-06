{
  description = "ClassyQ";
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
    flake-utils.url = "github:numtide/flake-utils";
    mach-nix.url = "github:DavHau/mach-nix/3.4.0";
  };

  outputs = { self, nixpkgs, flake-utils, mach-nix }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
        pythonEnv = mach-nix.lib."${system}".mkPython {
          requirements = builtins.readFile ./docs/requirements.txt +  ''
            jupyterlab
            sympy
            numpy
          '';
        };
      in
      {
        devShell = pkgs.mkShell.override { stdenv = pkgs.llvmPackages_14.stdenv; } {
          nativeBuildInputs = with pkgs; [
            clang-analyzer
            clang-tools
            cmake
            doctest
	    doxygen
            eigen
            fmt_8
            hdf5
            highfive
            lldb
            llvmPackages_14.clang-manpages
            llvmPackages_14.openmp
            ninja
            pythonEnv
            spdlog
            valgrind
          ];

          hardeningDisable = [ "all" ];

          NINJA_STATUS = "[Built edge %f of %t in %e sec] ";
        };
      });
}
