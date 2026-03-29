include(GNUInstallDirs)

# Installs .py under ${CMAKE_INSTALL_PREFIX}/python/... (layout relative to
# py_script_dir), and symlinks in BINDIR (stem name) to those files.
macro(smesh_install_python_scripts py_script_dir)
    if(IS_ABSOLUTE "${py_script_dir}")
        set(PY_SCRIPT_DIR "${py_script_dir}")
    else()
        set(PY_SCRIPT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${py_script_dir}")
    endif()

    file(GLOB_RECURSE PY_SCRIPTS CONFIGURE_DEPENDS
        "${PY_SCRIPT_DIR}/*.py"
    )

    foreach(py IN LISTS PY_SCRIPTS)
        file(RELATIVE_PATH rel_path "${PY_SCRIPT_DIR}" "${py}")
        get_filename_component(py_name "${py}" NAME)
        get_filename_component(py_stem "${py}" NAME_WE)
        get_filename_component(rel_dir "${rel_path}" DIRECTORY)

        if(rel_dir STREQUAL "")
            set(_py_dest "python")
        else()
            set(_py_dest "python/${rel_dir}")
        endif()

        install(PROGRAMS "${py}" DESTINATION "${_py_dest}")

        install(CODE
            "set(_root \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}\")
             set(_bindir \"\${_root}/${CMAKE_INSTALL_BINDIR}\")
             set(_pyfile \"\${_root}/${_py_dest}/${py_name}\")
             file(RELATIVE_PATH _rel \"\${_bindir}\" \"\${_pyfile}\")
             file(CREATE_LINK \"\${_rel}\" \"\${_bindir}/${py_stem}\" SYMBOLIC RESULT _res)
             if(NOT _res STREQUAL \"0\")
               message(FATAL_ERROR
                 \"Failed to create symlink: \${_bindir}/${py_stem} -> \${_rel}: \${_res}\")
             endif()"
        )
    endforeach()
endmacro()
