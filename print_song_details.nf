'''
This script takes two types of inputs: (a) Channels and (b) List and prints them using the workflow statement.
The stdout makes sure the output is printed in a terminal.
'''


ch1 = Channel.of('Satranga','Bekhayali','Pehli Dafa')
ch2 = Channel.of('Arijit Singh','Vishal Mishra','Atif Aslam')
list = ['\n This is a list of top 3 bollywood songs and their artists.', '(Reviews are subject to change)']

process song_artist{
    input:
    val x
    val y

    output:
    stdout

    script:
    """
    echo $x':' $y
    """
}

workflow {
    song_artist(ch1, ch2).view()
Channel
    .fromList(list)
    .view()
}

'''
A simplified way of writing this script would be to use a dictionary as input.
'''

song_details = [ "Satranga": "Arijit Singh", "Bekhayali":"Vishal Mishra", "Pehli Dafa": "Atif Aslam" ]
ch1 = Channel.of(song_details)

workflow {
    ch1.view()
}



'''
Modifying and adding the output to an external file will save the name of the top song as follows:
'''


ch1 = Channel.of('Satranga','Bekhayali','Pehli Dafa')
ch2 = Channel.of('Arijit Singh','Vishal Mishra','Atif Aslam')
list = ['\n This is a list of top 3 bollywood songs and their artists.', '(Reviews are subject to change)']

process song_artist{
    input:
    val x
    val y

    output:
    stdout

    script:
    """
    echo $x':' $y > all_songs
    cat all_songs | tail -n 1 | head -c 12 > best_song.txt

    """
}

workflow {
    song_artist(ch1, ch2).view()
Channel
    .fromList(list)
    .view()
}