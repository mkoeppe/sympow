SYMPOW Data setup
=================

__SYMPOW__, as forked, runs as any regular executable, that is to say,
not in an esoteric manner: intermediate computing steps that were meant
to be executed _by hand_ by end-users are now launched on the fly.
It also allows one to use a larger set of precomputed data which are
now managed system wide rather than locally.

The user-computed plain data mesh files and their binary counter parts,
are placed by default in the cache directory `HOME/.sympow` . This
default/ historical cache folder can be overridden through the environment
variable `SYMPOW_CACHEDIR`: the basename of the so passed path is assumed
to begin with `sympow`, otherwise the effective path is assumed to be
`SYMPOW_CACHEDIR/sympow`.

If the cache folder `SYMPOW_CACHEDIR` (or `SYMPOW_CACHEDIR/sympow`) does not
exists, then it is created on the fly with respect to permissions and
privileges whenever the following conditions are satisfied:
 (i) it terminates with at least three directory separators (the creation
   is group-centric when the number of separators is greater or equal to six,
   user-centric otherwise);
 (ii) its parent directory exists.
Along the same vain, the default/historical cache directory `HOME/.sympow` will
be created on the fly if it is nonexistent (the creation is user-centric).

The precomputed plain data mesh files are looked in `/usr/local/share/sympow`;
this default directory can be overridden through the environment variable
`SYMPOW_PKGDATADIR`. Their binary counter parts, the binary data mesh files,
are created on the fly and stored in the cache directory `/var/cache/sympow`;
this default directory can be overrident through the environment variable
SYMPOW_PKGCACHEDIR. The scripts effectivelly employed to compute data can be
found in `/usr/local/lib/sympow`; this default directory can be overridden through
the environment variable `SYMPOW_PKGLIBDIR`. Please note that the environment
variables introduced here (`SYMPOW_PKGDATADIR`, `SYMPOW_PKGLIBDIR`, and `SYMPOW_PKGCACHEDIR`)
are meant for advanced usage or debugging.

The data in `SYMPOW_CACHEDIR` and `SYMPOW_CACHEDIR/sympow` are authoritative.

For minutes details, you want to peruse the associated patches as imported
from Debian (package version 1.023-9) in release 2.023.1: _system wide scheme_
and _on fly new data_, commits f9add3 and c6f6cf, respectively.

Note that the `/usr/local` and `/var` prefixes employed above are the defaults
for the environment variables `PREFIX` and `VARPREFIX`, respectively, used by
Configure at configuration time.

