class Pfnet < Formula
  desc "Library for modeling and analyzing electric power networks"
  homepage "https://github.com/ttinoco/PFNET"
  head "https://github.com/ttinoco/PFNET.git"
  url "https://github.com/ttinoco/PFNET/archive/1.2.1.tar.gz"
  sha256 "ba0ddfb6eacb617d6a4ae7e425305b2a982d921f9b0cb706a7fb14a6a3214dd2"
  
  depends_on :python => :recommended if MacOS.version <= :snow_leopard
  depends_on :python3 => :recommended
  depends_on "numpy" => :python
  depends_on "graphviz" => :recommended

  def install
    args = ["NO_RAW_PARSER=1"]
    args << "NO_GRAPHVIZ=1" if build.without? "graphviz"

    pyargs = ["--no_raw_parser"]
    pyargs << "--no_graphviz" if build.without? "graphviz"

    system "make", *args


    Language::Python.each_python(build) do |python, version|
      dest_path = lib/"python#{version}/site-packages"
      dest_path.mkpath
      
      cd "python"

      system python, "setup.py", "build", *pyargs, "install", "--prefix=#{prefix}"

      cd buildpath
      lib.install "lib/libpfnet.so"
    end
  end

  def caveats
    if build.with?("python") && !Formula["python"].installed?
      homebrew_site_packages = Language::Python.homebrew_site_packages
      user_site_packages = Language::Python.user_site_packages "python"
      <<-EOS.undent
        If you use system python (that comes - depending on the OS X version -
        with older versions of numpy, scipy and matplotlib), you may need to
        ensure that the brewed packages come earlier in Python's sys.path with:
          mkdir -p #{user_site_packages}
          echo 'import sys; sys.path.insert(1, "#{homebrew_site_packages}")' >> #{user_site_packages}/homebrew.pth
      EOS
    end
  end

  test do
    Language::Python.each_python(build) do |python, _version|
      system python, "-c", "import pfnet; assert True"
    end
  end
end
