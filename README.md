## About the workflow of paper

### Basics of Git
we do not use git branches, we push everything into the master branch directly, so the following tutorial might be enough.

https://javascript.plainenglish.io/learn-git-in-less-than-10-minutes-62009d9416

Feel free to open an issue if you have any trouble using git features.

### To compile

1. Download this repo with `git clone`,
2. Open the main folder of this project,
3. Open a terminal and type
```bash
$ cd paper
$ latexmk -pdf paper.tex
```
4. Open the generated file `paper.pdf`

### To upload changes
```bash
$ git add file1 file2   # a lazy alternative is `git add -A` to add all changes
$ git status     # check status, make sure you add the wrong files.
```
NOTE: `# ...` is a comment, you do not need to type this part.

If you added the correct file, commit and upload changes
```bash
$ git commit -m 'some comments'
$ git push origin master
```
Otherwise, type
```bash
$ git reset
```
and restart from begining.

#### To write comments
1. open an issue in this repository
2. write comments directly in the file `paper/paper.tex`
    * \red{contents...}  # Shengtao Wang
    * \greed{contents...}  # Xun Gao

### To check the change history
If you want to check the difference with the version 2 commits backward, just type
```bash
$ cd paper
$ ./tdiff 2
```
If you see any warning or error, ignore and press enter.
To check the commit history, type

```bash
git log
```