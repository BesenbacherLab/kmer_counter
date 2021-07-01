# kmer_counter

[![ci](https://github.com/besenbacher/kmer_counter/workflows/ci/badge.svg)](https://github.com/besenbacher/kmer_counter/actions?query=workflow%3Aci)
[![documentation](https://img.shields.io/badge/docs-mkdocs%20material-blue.svg?style=flat)](https://besenbacher.github.io/kmer_counter/)
[![pypi version](https://img.shields.io/pypi/v/kmer_counter.svg)](https://pypi.org/project/kmer_counter/)
[![gitter](https://badges.gitter.im/join%20chat.svg)](https://gitter.im/kmer_counter/community)

Count kmers in regions or at SNVs or at indel breakpoints.

## Requirements

kmer_counter requires Python 3.6 or above.

<details>
<summary>To install Python 3.6, I recommend using <a href="https://github.com/pyenv/pyenv"><code>pyenv</code></a>.</summary>

```bash
# install pyenv
git clone https://github.com/pyenv/pyenv ~/.pyenv

# setup pyenv (you should also put these three lines in .bashrc or similar)
export PATH="${HOME}/.pyenv/bin:${PATH}"
export PYENV_ROOT="${HOME}/.pyenv"
eval "$(pyenv init -)"

# install Python 3.6
pyenv install 3.6.12

# make it available globally
pyenv global system 3.6.12
```
</details>

## Installation

With `pip`:
```bash
python3.6 -m pip install kmer_counter
```

With [`pipx`](https://github.com/pipxproject/pipx):
```bash
python3.6 -m pip install --user pipx

pipx install --python python3.6 kmer_counter
```
