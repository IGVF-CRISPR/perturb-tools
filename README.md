# single-cell-crispr-tool
Utility package for single-cell CRISPR screens pipeline with paired single-cell omics and guide counts

### Planning

We have a **[PR](https://github.com/IGVF-CRISPR/sccrispr-tool/pull/1)** for planning laying out modules and sub-modules for key functionalities required for the first-draft pipeline.

In general the structure could look something like the following:

```
IGVF-CRISPR/single-cell-crispr-pipeline/
│
├── LICENSE
├── notebooks
│   ├── examples
│   
├── README.md
├── requirements.txt
├── setup.py
│
├── sccrispr-tool/
│   ├── __init__.py
│   ├── data_structure/
│   │    ├── ...
│   │    
│   ├── qc/
│   │    ├── ...
│   │     
│   ├── _external_tools/
│   │    ├── .../
│   │    
```

Adhering to the above structure (or some structure - as long as it follows somewhat closely to the [PEP8 style guide](https://www.python.org/dev/peps/pep-0008/) will enable seamless collaboration regardless of the module you're working on - this is especially important for code review. 

### Contributing

In order to stay organzed, let's all contribute through PRs. Think of a PR as the main topic you are planning to contribute (`qc` or `guide_counting`). If PRs and issues are foreign to you, just ask! The best way to learn git workflows is through doing. Every time we open a PR, we should organize sub-tasks as issues and link them to that PR. Conversations, feedback, and requested changes can all be mediated through the PRs. 

#### Acknowledgement
`README`/collaborative package development practice advice from @mvinyard

