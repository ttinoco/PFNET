using BinDeps

@BinDeps.setup

pfnetname = "pfnet-1.3.3"
libpfnet = library_dependency("libpfnet")
prefix=joinpath(BinDeps.depsdir(libpfnet),"usr")
srcdir=joinpath(BinDeps.depsdir(libpfnet),"src")
downloads=joinpath(BinDeps.depsdir(libpfnet),"downloads")

# Linux
if is_linux()
    provides(SimpleBuild,
             (@build_steps begin
                ChangeDirectory(downloads)
                @build_steps begin
                  `gunzip -kf $pfnetname.tar.gz`
                  `tar -xvf $pfnetname.tar`
                  `rm $pfnetname.tar`
                  `mkdir -p $srcdir`
                  `rm -rf $srcdir/$pfnetname`
                  `mv $pfnetname $srcdir`
                end
              end,
              @build_steps begin
                ChangeDirectory(joinpath(srcdir, pfnetname))
                @build_steps begin
                  `./configure --prefix=$prefix`
                  `make clean`
                  `make uninstall`
                  `make`
                  `make check`
                  `make install`
                end
              end),
             libpfnet)
end

# Apple
if is_apple()
    provides(SimpleBuild,
             (@build_steps begin
                ChangeDirectory(downloads)
                @build_steps begin
                  `gunzip -kf $pfnetname.tar.gz`
                  `tar -xvf $pfnetname.tar`
                  `rm $pfnetname.tar`
                  `mkdir -p $srcdir`
                  `rm -rf $srcdir/$pfnetname`
                  `mv $pfnetname $srcdir`
                end
              end,
              @build_steps begin
                ChangeDirectory(joinpath(srcdir, pfnetname))
                @build_steps begin
                  `./configure --prefix=$prefix LDFLAGS="-Wl,-install_name,@loader_path/libpfnet.dylib"`
                  `make clean`
                  `make uninstall`
                  `make`
                  `make check`
                  `make install`
                end
              end),
             libpfnet)
end

# Windows
if is_windows()

end

@BinDeps.install Dict(:libpfnet => :libpfnet)
