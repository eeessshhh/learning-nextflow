'''
Using the shebang operator, you can define the path to the python tool located within the system.
Python commands can then be included in the script block of the process.
'''

process run_python_script {

    output:
    stdout

    script:
    """
    #!/bin/python3

    x = 'Teenage Mutant'
    y = 'Ninja Turtles'
    print("These are the outputs:", x,y)
    """
}

workflow {
    run_python_script()
}

