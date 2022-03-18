{
  description = "ClassyQ";
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
    flake-utils.url = "github:numtide/flake-utils";
    pypi-deps-db = {
      url = "github:DavHau/pypi-deps-db";
      flake = false;
    };
    mach-nix = {
      url = "mach-nix/3.4.0";
      inputs.nixpkgs.follows = "nixpkgs";
      inputs.flake-utils.follows = "flake-utils";
      inputs.pypi-deps-db.follows = "pypi-deps-db";
    };
  };

  outputs = { self, nixpkgs, flake-utils, mach-nix, pypi-deps-db }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
        pythonEnv = mach-nix.lib."${system}".mkPython {
          requirements = ''
            jupyterlab
            numpy
          '';
        };
      in
      {
        devShell = pkgs.mkShell.override { stdenv = pkgs.llvmPackages_13.stdenv; } {
          nativeBuildInputs = with pkgs; [
            catch2
            clang-analyzer
            clang-tools
            cmake
            doctest
            eigen
            fmt_8
            hdf5
            highfive
            llvmPackages_13.openmp
            ninja
            spdlog
          ];
          buildInputs = with pkgs; [
            pythonEnv
            pythonEnv.pkgs.jax
            pythonEnv.pkgs.jaxlib
            lldb
            valgrind
          ];

          hardeningDisable = [ "all" ];

          NINJA_STATUS = "[Built edge %f of %t in %e sec] ";
        };
      });
}
