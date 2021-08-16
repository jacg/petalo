test:
	cargo test
	cargo test -p cmlem
	pytest -v python bindings
