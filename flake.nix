{
  description = "ClassyQ";
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
      in
      {
        devShell = pkgs.mkShell.override { stdenv = pkgs.llvmPackages_14.stdenv; } {
          nativeBuildInputs = with pkgs; [
            clang-analyzer
            clang-tools
            cmake
            doctest
            eigen
            fmt_8
            hdf5
            highfive
            lldb
            llvmPackages_14.clang-manpages
            llvmPackages_14.openmp
            ninja
            spdlog
            valgrind
          ];

          hardeningDisable = [ "all" ];

          NINJA_STATUS = "[Built edge %f of %t in %e sec] ";
        };
      });
}
