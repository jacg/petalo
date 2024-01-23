# -*-Makefile-*-

test colours='':
	#!/usr/bin/env sh
	FAILED=0
	just test-rust   {{colours}} || FAILED=1
	just test-python {{colours}} || FAILED=2
	exit $FAILED

test-rust colours='':
	cargo nextest run {{colours}} --workspace

test-python colours='': python-build-bindings
	#!/usr/bin/env sh
	FAILED=0
	just test-python-pure     {{colours}} || FAILED=1
	just test-python-bindings {{colours}} || FAILED=2
	exit $FAILED

test-python-pure colours='':
	pytest -v {{colours}} src


test-python-bindings colours='': python-build-bindings
	pytest -v {{colours}} bindings



python-build-bindings profile='default':
	#!/usr/bin/env sh
	cd bindings
	case {{profile}} in
		default )
			cargo build
			ln -fs ../target/debug/libfulano.so fulano.so
			;;
		release )
			cargo build --release
			ln -fs ../target/release/libfulano.so fulano.so
			;;
		* )
			echo Unknown profile {{profile}}
			exit 1
	esac
