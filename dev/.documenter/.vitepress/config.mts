import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import mathjax3 from "markdown-it-mathjax3";
import footnote from "markdown-it-footnote";

// https://vitepress.dev/reference/site-config
export default defineConfig({
  base: '/BioMarkovChains.jl/',// TODO: replace this in makedocs!
  title: 'BioMarkovChains.jl',
  description: "A VitePress Site",
  lastUpdated: true,
  cleanUrls: true,
  outDir: '../final_site', // This is required for MarkdownVitepress to work correctly...
  
  ignoreDeadLinks: true,

  markdown: {
    math: true,
    config(md) {
      md.use(tabsMarkdownPlugin),
      md.use(mathjax3),
      md.use(footnote)
    },
    theme: {
      light: "github-light",
      dark: "github-dark"}
  },
  themeConfig: {
    outline: 'deep',
    
    search: {
      provider: 'local',
      options: {
        detailedView: true
      }
    },
    nav: [
{ text: 'Home', link: '/index' },
{ text: 'Towards Markov Chains', link: '/biomarkovchains' },
{ text: 'API', link: '/api' }
]
,
    sidebar: [
{ text: 'Home', link: '/index' },
{ text: 'Towards Markov Chains', link: '/biomarkovchains' },
{ text: 'API', link: '/api' }
]
,
    editLink: { pattern: "https://https://camilogarciabotero.github.io/BioMarkovChains.jl/edit/main/docs/src/:path" },
    socialLinks: [
      { icon: 'github', link: 'https://camilogarciabotero.github.io/BioMarkovChains.jl' }
    ],
    footer: {
      message: 'Made with <a href="https://luxdl.github.io/DocumenterVitepress.jl/dev/" target="_blank"><strong>DocumenterVitepress.jl</strong></a><br>',
      copyright: `© Copyright ${new Date().getUTCFullYear()}.`
    }
  }
})
