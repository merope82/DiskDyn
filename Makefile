.PHONY: all clean tgz

all:
	@make --no-print-directory -C src
	@make --no-print-directory -C tools

clean:
	@make --no-print-directory -C src clean
	@make --no-print-directory -C tools clean

nuke:
	@make --no-print-directory -C src nuke
	@make --no-print-directory -C tools nuke
	rm -f bin/*

tgz:
	tar cf diskdyn.tar ../DiskDyn_upload
	gzip diskdyn.tar
