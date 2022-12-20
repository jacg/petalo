# -*-Makefile-*-

test: test-rust test-python hand

hand: ms mc ss sc

mc:
	cargo run -q --bin make_simset_lors -- ~/src/simpet/samples/NEMA_IEC_scatter_counts/custom_det_hf.hist   -t custom    -n 10

ms:
	cargo run -q --bin make_simset_lors -- ~/src/simpet/samples/NEMA_IEC_scatter_counts/standard_phg_hf.hist -t standard  -n 10

ss:
	cargo run -q -p         simset      -- ~/src/simpet/samples/NEMA_IEC_scatter_counts/standard_phg_hf.hist -t standard  -n 5

sc:
	cargo run -q -p         simset      -- ~/src/simpet/samples/NEMA_IEC_scatter_counts/custom_det_hf.hist   -t custom    -n 10

test-rust colours='':
	just test-rust-pure {{colours}}

test-rust-pure colours='':
	cargo nextest run {{colours}} --workspace --exclude bindings


test-python colours='': python-build-bindings
	pytest -v {{colours}} src bindings


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
