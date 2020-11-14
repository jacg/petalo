# ----- Install home-manager ------------------------------------

# Default location where home-manager configuration will be placed. This can be
# overridden by calling this script thus:
#
# HM_DIR=my-choice-of-home-manager-location sh install-home-manager-and-maybe-spacemacs.sh
: ${HM_DIR:=$HOME/my-home-manager}

echo home-manager configuration will be placed in $HM_DIR


# Download and unpack the home-manager config template
cd /tmp
# TODO: change URL once this branch is merged!
curl -L https://github.com/jacg/petalo/tarball/devrust-petalo > /tmp/petalo.tgz
#curl -L https://github.com/jacg/petalo/tarball/master > /tmp/petalo.tgz
nix-shell -p pkgs.gnutar --run "tar xvf /tmp/petalo.tgz --wildcards '*/nix/home-manager' --strip-components=2"
mv home-manager $HM_DIR

# Bootstrap your personal home-manager installation and configuration
cd $HM_DIR
nix-shell bootstrap-home-manager
# You should also run `nix-shell bootstrap-home-manager` (from $HM_DIR) whenever
# you change the home-manager version in ./sources.nix

# ----- Simple check of home-manager installation ------------------
echo "Checking home-manager's ability to install new packages"

# First confirm that `cowsay` is not installed:
cowsay "Is cowsay installed?" || echo "Nope, cowsay is not installed yet"

# Now edit your `home.nix` to include pkgs.cowsay in `home.packages`
nix-shell -p pkgs.gnused --run "sed -i 's/#pkgs.cowsay/pkgs.cowsay/g' $HM_DIR/nixpkgs/home.nix"

# Instruct home-manager to switch to your new configuration
home-manager -b bck switch

# And check that cowsay is installed now
cowsay "Home manager has managed to install cowsay"

# ------ Install space macs with Rust development tools -----

# By default, spacemacs will be installed with emacs editing style as the
# default
: ${SPACEMACS:=emacs}

# You can choose not to install spacemacs at all:
#
#    SPACEMACS=no sh install-home-manager-and-maybe-spacemacs.sh

# You can choose to set the vim editing style as default:
#
#    SPACEMACS=vim sh install-home-manager-and-maybe-spacemacs.sh

# TODO: provide option to install spacemacs ALONGSIDE "normal" emacs

case "$SPACEMACS" in
    no)
        echo "NOT installing Spacemacs with Rust devtools"
        exit 0
        ;;
    emacs|vim)
        echo "Installing Spacemacs with Rust devtools"
        echo "Spacemacs default editing style: $SPACEMACS"
        ;;
    *)
        echo "ERROR: invalid SPACEMACS value: $SPACEMACS"
        echo "Allowed values: emacs, vim, no"
        exit 1
        ;;
esac

# Uncomment the relevant parts of home.nix
nix-shell -p pkgs.gnused --run "sed -i 's/# <uncomment-to-install-spacemacs> # //g' $HM_DIR/nixpkgs/home.nix"

# Instruct home-manager to install the new configuration (including spacemacs)
home-manager -b bck switch

# Get spacemacs to install emacs packages
emacs -batch --eval="(setf (symbol-function 'y-or-n-p) (lambda (p) (print p) t))" -l ~/.emacs.d/init.el

# TODO maybe clear .emacs.d/{.cache,elpa} before starting?
