function flag_make_doc(path)

cd(path)
m2html('mfiles', 'src/main/matlab', 'htmldir', 'doc/matlab');

end