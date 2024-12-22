{
  description = "Wave equation Solver";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = nixpkgs.legacyPackages.${system};

        pythonEnv = pkgs.python3.withPackages
          (ps: with ps; [ numpy matplotlib sympy tqdm scipy ]);

        cppProject = pkgs.stdenv.mkDerivation {
          pname = "wave_solver";
          version = "0.1.0";
          src = ./cpp;

          nativeBuildInputs = with pkgs; [
            cmake
            ninja
            pkg-config
          ];

          buildInputs = with pkgs; [
            gnuplot
          ];

          cmakeFlags = [
            "-DCMAKE_BUILD_TYPE=Release"
            "-GNinja"
          ];

          installPhase = ''
            mkdir -p $out/bin
            cp wave_demo $out/bin/
          '';
        };

        pythonProject = pkgs.stdenv.mkDerivation {
          pname = "python-boundary_problem";
          version = "0.1.0";
          src = ./py;

          nativeBuildInputs = [ pythonEnv ];

          installPhase = ''
            mkdir -p $out/bin $out/lib/python
            cp -r . $out/lib/python/
            cat > $out/bin/wave_solver_py <<EOF
            #!${pythonEnv}/bin/python
            import sys
            sys.path.insert(0, '$out/lib/python')
            from demo import main
            main()
            EOF
            chmod +x $out/bin/wave_solver_py
          '';
        };

      in {
        packages = {
          cpp = cppProject;
          py = pythonProject;
          default = cppProject;
        };

        apps = {
          cpp = flake-utils.lib.mkApp {
            drv = cppProject;
            name = "wave_demo";
          };
          py = flake-utils.lib.mkApp {
            drv = pythonProject;
            name = "wave_solver_py";
          };
          default = self.apps.${system}.cpp;
        };

        devShells.default = pkgs.mkShell {
          packages = with pkgs; [
            # Build tools
            cmake
            ninja
            gnumake
            pkg-config
            ccache

            # Development tools
            gdb
            valgrind
            bear
            clang-tools

            # C++ standard library and compiler
            gcc
            stdenv.cc.cc.lib

            # Libraries
            boost
            gnuplot

            # Python environment
            pythonEnv
            pyright
          ];

          shellHook = ''
            export CCACHE_DIR=$HOME/.ccache
            export PATH="$HOME/.ccache/bin:$PATH"

            # Add C++ standard library headers to CPATH
            export CPATH="${pkgs.gcc-unwrapped}/include/c++/${pkgs.gcc-unwrapped.version}:${pkgs.gcc-unwrapped}/include/c++/${pkgs.gcc-unwrapped.version}/${pkgs.stdenv.targetPlatform.config}:$CPATH"

            alias c=clear

            echo "======================================"
            echo "$(cmake    --version | head -n 1)"
            echo "$(ninja    --version | head -n 1)"
            echo "$(g++      --version | head -n 1)"
            echo "$(python   --version | head -n 1)"
            echo "$(gnuplot  --version | head -n 1)"
            echo ""
            echo "Build the project:  nix build"
            echo "Run C++ project:    nix run .#cpp"
            echo "Run Python project: nix run .#py"
            echo ""
            echo "Happy coding!"
          '';
        };
      });
}
