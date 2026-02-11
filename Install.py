from tkinter import *

titlef = ["Malgun Gothic", 17, "bold"]
textf = ["Malgun Gothic", 15]
win = Tk()
win.title("How to Install")
l1 = Label(win, text="How to install Universe Explorer?", font=titlef)
l1.pack()

l2 = Label(win, text="(ENG) Universe Explorer requires Python and several Python libraries.\n" \
"When installing Python, be sure to check the Add python.exe to PATH checkbox.\n" \
"Enter the following commands in order in the terminal window.\n" \
"If the terminal outputs that the library is already installed, just skip it.", font=textf)
l2.pack()

l3 = Label(win, text="(KOR) Universe Explorer를 실행하려면 Python과 여러 Python 라이브러리들이 필요합니다.\n" \
"Python을 설치할 때 Add python.exe to PATH 체크박스를 꼭 체크하고 설치해 주세요.\n" \
"아래 명령어들을 터미널 창에 차례대로 입력 해 주세요.\n" \
"터미널에서 라이브러리가 이미 설치되어 있다고 출력되면 그냥 넘어가세요.", font=textf)
l3.pack()

t1 = Text(win, width=50, height=10, font=textf)
t1.insert("1.0", "pip install sv_ttk\npip install pillow\n" \
"pip install requests\npip install astropy\n" \
"pip install numpy\npip install matplotlib\n" \
"pip install astroquery\npip install lightkurve")
t1.pack()

win.mainloop()
