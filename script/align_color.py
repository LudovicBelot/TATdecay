import sys

def main():

    file2color = sys.argv[1]
    
    str_html_ali = ali2color(file2color)


def ali2color(file):

    n_index = 0
    str_html = "<!DOCTYPE html>\n<html>\n<body>"


    with open(file, "r") as f:
        for line in f:
            n_index += 1
            n_index_aa = 0
            if n_index == 1:
              str_html += f"<b>{line}</b>\n"
            elif n_index == 2:
                str_html += "<p style='line-height:0.3em;'>"+line.split('\t\t')[0].strip() + "&emsp;&emsp;"+"<span style='background-color: green'>" + line.split('\t\t')[1].strip() +"</span></p>\n"
                ref_seq = line.split('\t\t')[1]
            else :
                if line.startswith(" ") == False:
                    str_html += "<p style='line-height:0.3em;'>"+ line.split('\t\t')[0].strip() + "&emsp;&emsp;"
                    genome_name = line.split('\t\t')[0].strip()
                else :
                    str_html += "<p style='line-height:0.3em;'>"+ "<span style='color: white'>" + genome_name +"</span>" + "&emsp;&emsp;"

                for aa in line.split('\t\t')[1]:
                    if aa != "\n":
                        str_html += color_aa(aa, ref_seq[n_index_aa])
                        n_index_aa += 1
                    else :
                        str_html += "</span></p>\n"

    str_html += "</html>\n</body>"

    with open(f"{file.rsplit('.',1)[0]}.html", "w") as f:
        f.write(str_html)


def color_aa(new_aa, ref_aa):
    
    if new_aa == " ":
        return "<span style='color: white'>" + ref_aa +"</span>"
    else :
        if new_aa == ref_aa:
            return "<span style='background-color: green'>" + new_aa +"</span>"
        else :
            if new_aa == "*" or new_aa == "-":
                return "<span style='background-color: red'>" + new_aa +"</span>"
            else :
                return "<span style='background-color: grey'>" + new_aa +"</span>"


#<span style="color: brown">WORD</span>



if __name__ == "__main__":
    main()