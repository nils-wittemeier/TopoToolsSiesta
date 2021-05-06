#!/bin/bash
function install_test_sisl {
    echo ""
    echo " Will try and run sisl"
    echo "    import sisl ; print(sisl.geom.graphene())"
    python3 -c "import sisl ; print(sisl.geom.graphene())"
    if [ $? -ne 0 ]; then
        echo "Failed running sisl, please mail the organizer with the error message (unless some of the installations failed)"
    fi
}

function install_test_z2pack {
    echo ""
    echo " Will try and z2pack after patch"
    echo "    import z2pack; z2pack.hm.System(hamilton=(lambda k: [[0]]), basis_overlap=(lambda k: [[1]]))"
    python3 -c "import z2pack; z2pack.hm.System(hamilton=(lambda k: [[0]]), basis_overlap=(lambda k: [[1]]))"
    if [ $? -ne 0 ]; then
        echo "Failed running Z2Pack, please mail the organizer with the error message (unless some of the installations failed)"
    fi
}

install_test_sisl
install_test_z2pack
