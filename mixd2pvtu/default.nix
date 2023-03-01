{
    pkgs ? import (builtins.fetchTarball {

      # Descriptive name to make the store path easier to identify
      name = "nixpkgs-unstable-2022-03-17";

      # Commit hash for nixos-unstable as of 2018-09-12
      url = "https://github.com/nixos/nixpkgs/archive/3eb07eeafb52bcbf02ce800f032f18d666a9498d.tar.gz";

      # Hash obtained using `nix-prefetch-url --unpack <url>`
      sha256 = "1ah1fvll0z3w5ykzc6pabqr7mpbnbl1i3vhmns6k67a4y7w0ihrr";
    
    }) {}

}:

let

    version = "0.1";

    vtk910-mpi = (pkgs.vtk_9.overrideAttrs(oldAttrs: rec{

            majorVersion = "9.1";
            minorVersion = "0";

            version = "${majorVersion}.${minorVersion}";

            src = pkgs.fetchurl {
            url = "https://www.vtk.org/files/release/${majorVersion}/VTK-${version}.tar.gz";
            sha256 = "8fed42f4f8f1eb8083107b68eaa9ad71da07110161a3116ad807f43e5ca5ce96";
            };

            cmakeFlags = oldAttrs.cmakeFlags ++ [
                "-DVTK_USE_MPI=ON"
                "-DVTK_SMP_IMPLEMENTATION_TYPE=TBB"
            ];

            buildInputs = oldAttrs.buildInputs ++ [
                pkgs.tbb
                pkgs.mpi
            ];

            # The default script currently does something stupid and the
            # install breaks
            preFixup = "";

            postInstall = ''
            ln -s $out/include/vtk-${majorVersion} $out/include/vtk
            '';

    }) ).override { stdenv = pkgs.gcc11Stdenv; };

in pkgs.gcc11Stdenv.mkDerivation rec {
    pname = "mixd2pvtu";
    inherit version;

    src = ./src;

    nativeBuildInputs = with pkgs;[
        cmake
        mpi
    ];

    # NOTE: Why can't I just bundle all the rest as dependencies of VTK? 
    buildInputs = with pkgs; [
        vtk910-mpi
        mesa
        libGL
        libGLU
        xorg.libX11
        tbb
    ];

    depsTargetTarget = with pkgs;[
        vtk910-mpi
    ];

}
